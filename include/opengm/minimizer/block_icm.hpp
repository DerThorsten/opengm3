#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/from_gm_factory.hpp"

#include "opengm/minimizer/factor_icm.hpp"
#include "opengm/minimizer/icm.hpp"
#include "opengm/minimizer/brute_force_naive.hpp"

namespace opengm{

template<class GM>
class BlockIcm;

template<class GM>
class BlockGeneratorBase
{
public:
    using gm_type = GM;
    using block_icm_type = BlockIcm<GM>;
    using variables_set_type = std::set<std::size_t>;

    virtual bool generate(block_icm_type * block_icm) = 0;

};

template<class gm_type>
using BlockGeneratorFactoryBase = detail::FromGmFactoryBase<BlockGeneratorBase<gm_type>, gm_type>;


template<class block_generator_type>
using BlockGeneratorFactory = detail::FromGmFactory<
    block_generator_type,
    BlockGeneratorBase<typename block_generator_type::gm_type>,
    typename block_generator_type::gm_type,
    typename block_generator_type::settings_type
>;



template<class GM>
class SingleVarBlockGenerator : public BlockGeneratorBase<GM>
{
public:
    using gm_type = GM;
    using block_icm_type = BlockIcm<GM>;
    using variables_set_type = std::set<std::size_t>;

    struct settings_type {
    };

    SingleVarBlockGenerator(const gm_type & gm, const settings_type settings = settings_type())
    :   m_gm(gm),
        m_settings(settings)
    {

    }

    bool generate(block_icm_type * block_icm) override
    {
        std::size_t some_var = 0;
        auto improvment  = block_icm->solve_submodel(&some_var, &some_var+1);
        return false;
    }
private:
    const gm_type & m_gm;
    settings_type m_settings;
};






template<class GM>
class BlockIcm : public MinimizerCrtpBase<GM, BlockIcm<GM> >{
public:
    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;

    using conditioned_submodel_builder_type = detail::ConditionedSubmodelBuilder<gm_type>;
    using sub_gm_type = typename conditioned_submodel_builder_type::sub_gm_type; 



    using default_sub_minimizer_type = FactorIcm<sub_gm_type>;
    using default_sub_minimier_factory_type = MinimizerFactory<default_sub_minimizer_type>;

    using sub_minimizer_factory_base_type = MinimizerFactoryBase<sub_gm_type>;
    using sub_minimizer_factory_ptr_type = std::shared_ptr<sub_minimizer_factory_base_type>;


    using block_generator_factory_base_type = BlockGeneratorFactoryBase<gm_type>;
    using block_generator_factory_ptr_type = std::shared_ptr<block_generator_factory_base_type>;

    using base_type::minimize;


    struct settings_type : public SolverSettingsBase{
        public:
            sub_minimizer_factory_ptr_type minimizer_factory;
            block_generator_factory_ptr_type block_generator_factory;
    };


    BlockIcm(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_settings(settings),
        m_labels(gm.num_variables(),0),
        m_energy(m_gm.evaluate(m_labels)),
        m_builder(m_gm)
    {
        if(! m_settings.block_generator_factory )
        {
            m_settings.block_generator_factory = std::make_shared<BlockGeneratorFactory<SingleVarBlockGenerator<gm_type>>>();
        }
        if(! m_settings.minimizer_factory )
        {
            m_settings.minimizer_factory = std::make_shared<MinimizerFactory<BruteForceNaive<sub_gm_type>>>();
        }
    }

    std::string name() const override{
        return "BlockIcm";
    }

    const gm_type & gm() const override{
        return m_gm;
    }

    const labels_vector_type & best_labels()override{
        return m_labels;
    }
    const labels_vector_type & current_labels() override{
        return m_labels;
    }

    value_type best_energy()  override{
        return m_energy;
    }
    value_type current_energy() override {
        return m_energy;
    }

    bool can_start_from_starting_point() override{
        return true;
    }
    void set_starting_point(const label_type  * labels)override{
        m_labels.assign(labels, labels + m_gm.num_variables());
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override{
        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);

        auto block_gen = m_settings.block_generator_factory->create(m_gm);
        while(true){
            auto continue_gen = block_gen->generate(this);
            // call the callback if solving the the submodel
            // improved the overall energy
            if(m_improvment)
            {
                callback();
            }

            if(!continue_gen)
            {
                break;
            }
        }
        m_energy = m_gm.evaluate(m_labels);

    }

    template<class VI_ITER>
    auto solve_submodel(VI_ITER free_vi_begin, VI_ITER free_vi_end)
    {
        m_improvment = false;
        auto delta = value_type(0);
        m_builder.condition(free_vi_begin, free_vi_end, m_labels, [&](auto && sub_gm, auto && current_sub_gm_labels){

            // solver the sub_gm
            std::cout<<"sub_gm "<<sub_gm.num_variables()<<" "<<sub_gm.num_factors()<<"\n";
            std::cout<<"create sub_gm_minimizer\n";
            auto sub_gm_minimizer = m_settings.minimizer_factory->create(sub_gm);
            std::cout<<"minimize\n";
            sub_gm_minimizer->minimize();

            std::cout<<"get best\n";
            const auto & sub_gm_labels = sub_gm_minimizer->best_labels();
            const auto & sub_gm_energy = sub_gm_minimizer->best_energy();
            const auto old_sub_energy  = sub_gm.evaluate(current_sub_gm_labels, current_sub_gm_labels+sub_gm.num_variables());
            delta = old_sub_energy - sub_gm_energy;
            if(delta > value_type(0))
            {
                m_improvment = true;

                for(auto sub_gm_vi=0; sub_gm_vi<sub_gm.num_variables(); ++sub_gm_vi)
                {
                    m_labels[m_builder.sub_gm_to_gm()[sub_gm_vi]] = sub_gm_labels[sub_gm_vi];
                }
            }

        });
        if(m_improvment)
        {
            m_energy -= delta;
        }
        return m_improvment;
    }
private:


    const GM & m_gm;
    settings_type m_settings;
    std::vector<label_type> m_labels;
    value_type m_energy;
    conditioned_submodel_builder_type m_builder;
    bool m_improvment;
};

template<class GM>
using BlockIcmFactory = MinimizerFactory<BlockIcm<GM>>;
}