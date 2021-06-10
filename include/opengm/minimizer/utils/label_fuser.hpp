#pragma once

#include <vector>
#include <memory>

#include "opengm/space.hpp"
#include "opengm/minimizer/factor_icm.hpp"
#include "opengm/graphical_model.hpp"

namespace opengm{

    template<class GM>
    class LabelFuser{
    public:

        using gm_type = GM;
        using label_type = typename gm_type::label_type;
        using value_type = typename gm_type::value_type;
        using label_vector_type = std::vector<label_type>;

        // the fm model
        using fuse_space = StaticNumLabelsSpace<label_type,2>;
        using fuse_gm_type = GraphicalModel<fuse_space, value_type>;
        using fuse_tensor_type = StaticNumLabelTensor<value_type, 2>;

        using default_fuse_minimizer_type = FactorIcm<fuse_gm_type>;
        using default_fuse_minimier_factory_type = MinimizerFactory<default_fuse_minimizer_type>;

        using fuse_minimizer_factory_base_type = MinimizerFactoryBase<fuse_gm_type>;
        using fuse_minimizer_factory_ptr_type = std::shared_ptr<fuse_minimizer_factory_base_type>;


        class settings_type
        {
        public:
            fuse_minimizer_factory_ptr_type minimizer_factory;
            value_type eps{1e-7};
        };

        LabelFuser(const gm_type & gm, const settings_type & settings = settings_type())
        :   m_gm(gm),
            m_settings(settings),
            m_gm_to_fuse_gm(m_gm.num_variables()),
            m_fuse_gm_to_gm(m_gm.num_variables()),
            m_in_fuse_gm(m_gm.num_variables()),
            m_fuse_gm_num_var(0),
            m_factor_labels(),
            m_fuse_factor_labels(),
            m_fuse_factor_vis(),
            m_free_var_pos(),
            m_fuse_gm(),
            m_fuse_gm_starting_point(m_gm.num_variables()),
            m_fuse_gm_unaries(m_gm.num_variables())
        {
            if(!m_settings.minimizer_factory)
            {
                m_settings.minimizer_factory = std::make_shared<default_fuse_minimier_factory_type>();
            }

            const auto max_arity = m_gm.max_arity();
            m_factor_labels.resize(max_arity);
            m_fuse_factor_labels.resize(max_arity);
            m_fuse_factor_vis.resize(max_arity);
            m_free_var_pos.resize(max_arity);
        }

        auto fuse(
            const label_vector_type & labels_a,
            label_vector_type & labels_b
        ){
            return this->fuse(labels_a, labels_b, labels_b);
        }

        auto fuse(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b,
            label_vector_type & fused
        ){
            // count num var and build mapping from / to fuse gm
            m_fuse_gm_num_var = setup_mapping(labels_a, labels_b);
            const auto energy_a = m_gm.evaluate(labels_a);
            if(m_fuse_gm_num_var > 0)
            {
                const auto energy_b = m_gm.evaluate(labels_b);

                // clear model and resize space
                m_fuse_gm.clear();
                m_fuse_gm.space().resize(m_fuse_gm_num_var);
                std::fill(m_fuse_gm_unaries.begin(), m_fuse_gm_unaries.begin() + m_fuse_gm_num_var, 0);

                for(auto && factor : m_gm)
                {
                    // setup variables of fuse factor
                    const auto fuse_factor_arity = setup_fuse_factor_variables(factor);

                    // fully included
                    if(fuse_factor_arity == factor.arity())
                    {
                        if(fuse_factor_arity == 1)
                        {
                            this->add_fully_included_unary_factor(labels_a, labels_b, factor);
                        }
                        else
                        {
                            this->add_fully_included_factor(labels_a, labels_b, factor);
                        }
                    }
                    // partially included
                    else if(fuse_factor_arity > 0)
                    {
                        if(fuse_factor_arity == 1)
                        {
                            this->add_unary_from_partially_included_factor(labels_a, labels_b, factor);
                        }
                        else
                        {
                            this->add_partially_included_factor(labels_a, labels_b, factor, fuse_factor_arity);
                        }
                    }
                }

                this->add_collected_unaries();

                // solve
                auto fuse_minimizer = m_settings.minimizer_factory->create(m_fuse_gm);
                if(fuse_minimizer->can_start_from_starting_point())
                {
                    std::fill(m_fuse_gm_starting_point.begin(), m_fuse_gm_starting_point.begin() + m_fuse_gm.num_variables(), energy_a < energy_b ? 0 :1);
                    fuse_minimizer->set_starting_point(m_fuse_gm_starting_point.data());
                }
                fuse_minimizer->minimize();
                const auto & fuse_gm_labels = fuse_minimizer->best_labels();

                // project the fuse-gm solution to the gm
                for(auto fvi=0; fvi<m_fuse_gm.num_variables(); ++fvi)
                {
                    const auto vi = m_fuse_gm_to_gm[fvi];
                    const auto fl =  fuse_gm_labels[fvi];
                    fused[vi] = (fl == 0 ? labels_a[vi] : labels_b[vi]);
                }

                // enforce that energy is not worse than any of candidates
                const auto e_fused = m_gm.evaluate(fused);
                if(e_fused >  std::min(energy_a, energy_b))
                {
                    fused = energy_a < energy_b ? labels_a : labels_b;
                    return std::min(energy_a, energy_b);
                }
                return e_fused;
            }
            else
            {
                fused = labels_a;
                return energy_a;
            }

        }

    private:

        void add_collected_unaries()
        {
            // add the unaries we collected
            for(auto fvi=0; fvi<m_fuse_gm_num_var; ++fvi)
            {
                const auto e0 = m_fuse_gm_unaries[fvi];
                if(std::fabs(e0) > m_settings.eps)
                {
                    auto tensor = std::make_unique<OptimizedBinaryUnary<value_type>>(e0);
                    m_fuse_gm.add_unary_factor(std::move(tensor), fvi);
                }
            }
        }

        template<class FACTOR>
        auto setup_fuse_factor_variables(FACTOR && factor){

            // check if the factors is  fully/partially/not included
            auto && vars = factor.variables();
            auto fuse_factor_arity = 0;
            for(auto ai=0; ai<factor.arity(); ++ai)
            {
                if(m_in_fuse_gm[vars[ai]])
                {
                    m_fuse_factor_vis[fuse_factor_arity] = m_gm_to_fuse_gm[vars[ai]];
                    m_free_var_pos[fuse_factor_arity] = ai;
                    ++fuse_factor_arity;
                }
            }
            return fuse_factor_arity;
        }


        template<class FACTOR>
        void add_fully_included_factor(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b,
            FACTOR && factor
        )
        {

            const auto arity = factor.arity();

            auto tensor = std::make_unique<fuse_tensor_type>(arity);

            auto && vars = factor.variables();
            detail::for_each_state<2>(arity, m_fuse_factor_labels, [&](auto && lables){

                // map from fuse-gm labels to gm labels
                for(std::size_t ai=0; ai<arity; ++ai)
                {
                    m_factor_labels[ai] = m_fuse_factor_labels[ai] == 0 ? labels_a[vars[ai]] : labels_a[vars[ai]];
                }

                // evaluate and assign factor values to tensor
                tensor->operator[](m_fuse_factor_labels.data()) = factor[m_factor_labels.data()];

            });

            // add the factor
            m_fuse_gm.add_factor(std::move(tensor), m_fuse_factor_vis.begin(), m_fuse_factor_vis.begin() + arity);
        }


        template<class FACTOR>
        void add_fully_included_unary_factor(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b,
            FACTOR && factor
        ){
            const auto vi = factor.variables()[0];
            const auto fvis = m_gm_to_fuse_gm[vi];
            const auto e0 = factor(labels_a[vi]);
            const auto e1 = factor(labels_b[vi]);
            m_fuse_gm_unaries[fvis] += (e0 - e1);
        }

        template<class FACTOR>
        void add_partially_included_factor(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b,
            FACTOR && factor,
            const std::size_t fuse_factor_arity
        )
        {
            const auto arity = factor.arity();
            auto && vars = factor.variables();


            auto tensor = std::make_unique<fuse_tensor_type>(fuse_factor_arity);

            for(std::size_t ai=0; ai<arity; ++ai)
            {
                // we can take the labels from labels_a or labels_b
                // here since m_factor_labels is only used
                // at placed i where labels_a[i] == labels_b[i]
                m_factor_labels[ai] = labels_a[vars[ai]];
            }

            detail::for_each_state<2>(fuse_factor_arity, m_fuse_factor_labels, [&](auto && labels){

                // translate from fuse factor labels to factor labels
                for(auto fai=0; fai<fuse_factor_arity; ++fai)
                {
                    // index wrt gm factor
                    const auto ai = m_free_var_pos[fai];
                    auto vi = vars[ai];

                    m_factor_labels[ai] = m_fuse_factor_labels[fai] == 0 ? labels_a[vi] : labels_b[vi];

                }
                tensor->operator[](m_fuse_factor_labels.data()) = factor[m_factor_labels.data()];
            });
            // add the factor
            m_fuse_gm.add_factor(std::move(tensor), m_fuse_factor_vis.begin(), m_fuse_factor_vis.begin() + fuse_factor_arity);
        }


        template<class FACTOR>
        void add_unary_from_partially_included_factor(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b,
            FACTOR && factor
        ){
            auto && vars = factor.variables();
            const auto pos = m_free_var_pos[0];
            const auto vi = vars[pos];
            const auto fvis = m_gm_to_fuse_gm[vi];

            for(auto ai=0; ai<factor.arity(); ++ai)
            {
                m_factor_labels[ai] = labels_a[vars[ai]];
            }

            const auto e0 = factor[m_factor_labels.data()];
            m_factor_labels[pos] = labels_b[vi];
            const auto e1 = factor[m_factor_labels.data()];

            m_fuse_gm_unaries[fvis] += (e0 - e1);
        }

        auto setup_mapping(
            const label_vector_type & labels_a,
            const label_vector_type & labels_b
        ){
            // count num var
            // and build mapping from / to fuse gm
            auto fuse_gm_num_var = 0;
            for(std::size_t vi=0; vi<m_gm.num_variables(); ++vi)
            {
                if(labels_a[vi] != labels_b[vi])
                {
                    m_gm_to_fuse_gm[vi] = fuse_gm_num_var;
                    m_fuse_gm_to_gm[fuse_gm_num_var] = vi;
                    m_in_fuse_gm[vi] = true;
                    ++fuse_gm_num_var;
                }
                else
                {
                    m_in_fuse_gm[vi] = false;
                }
            }
            return fuse_gm_num_var;
        }


        const GM & m_gm;
        settings_type m_settings;
        std::vector<std::size_t> m_gm_to_fuse_gm;
        std::vector<std::size_t> m_fuse_gm_to_gm;
        std::vector<bool> m_in_fuse_gm;
        std::size_t m_fuse_gm_num_var;

        std::vector<label_type> m_factor_labels;
        std::vector<label_type> m_fuse_factor_labels;
        std::vector<std::size_t> m_fuse_factor_vis;
        std::vector<std::size_t> m_free_var_pos;

        fuse_gm_type m_fuse_gm;
        std::vector<typename fuse_gm_type::label_type> m_fuse_gm_starting_point;
        std::vector<value_type> m_fuse_gm_unaries;
    };


}