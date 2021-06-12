#pragma once

#include <queue>

#include "opengm/opengm_config.hpp"
#include "opengm/datastructures/partition.hpp"
#include "opengm/graphical_model.hpp"
#include "opengm/factors_of_variables.hpp"

#include <gsl-lite/gsl-lite.hpp>




namespace opengm
{




namespace detail
{

    template<class GM>
    class TreeFinder
    {
    public:
        using gm_type = GM;
        using factors_of_variables_type = FactorsOfVariables<gm_type>;

    public:
        TreeFinder(const gm_type & gm, const factors_of_variables_type & factors_of_variables)
        :   m_gm(gm),
            m_factors_of_variables(factors_of_variables),
            m_neighbours(m_gm.num_variables())
        {
            // build simple nh graph
            for(auto vi=0; vi<m_gm.num_variables(); ++vi)
            {
                for(auto fi : m_factors_of_variables[vi])
                {
                    auto && factor = m_gm[fi];
                    if(factor.arity() > 1)
                    {
                        for(auto other_vi : factor.variables())
                        {
                            if(other_vi != vi)
                            {
                                m_neighbours[vi].insert(other_vi);
                                m_neighbours[other_vi].insert(vi);
                            }
                        }
                    }
                }
            }
        }


        void find(const std::size_t seed_var){

            opengm::Partition<std::size_t> ufd(m_gm.num_variables());

            std::set<std::size_t>   tree_vars;
            std::queue<std::size_t> queue;
            queue.push(seed_var);


            auto can_be_added = [&](auto candidate_vi){
                if(tree_vars.size() <= 1)
                {
                    return true;
                }
                else{
                    auto n_in_tree = 0;
                    for(auto neighbour_vi : m_neighbours[candidate_vi])
                    {
                        if(tree_vars.find(neighbour_vi) != tree_vars.end())
                        {
                            ++n_in_tree;
                            if(n_in_tree >= 2)
                            {
                                return false;
                            }
                        }
                    }
                    return true;
                }
            };


            while(!queue.empty())
            {
                // get candidate
                auto candidate_vi = queue.front();
                queue.pop();

                // only consider variables which are not
                // yet in the tree
                if(tree_vars.find(candidate_vi) == tree_vars.end())
                {
                    // check if the variable can be added
                    // without adding a cycle
                    if(can_be_added(candidate_vi))
                    {
                        tree_vars.insert(candidate_vi);
                        for(auto neighbour_vi : m_neighbours[candidate_vi])
                        {
                            queue.push(neighbour_vi);
                        }
                    }
                }
            }
            //std::cout<<tree_vars.size()<<"\n";
        }
    private:
        const gm_type & m_gm;
        const factors_of_variables_type & m_factors_of_variables;

        std::vector<std::set<std::size_t>> m_neighbours;
    };




    
    template<class GM>
    class ConditionedSubmodelBuilder
    {
    public:
        using gm_type = GM;
        using value_type = typename gm_type::value_type;
        using labels_vector_type = typename gm_type::labels_vector_type;


        // submodel
        using subspace_type = typename gm_type::space_type::subspace_type;
        using sub_gm_type = GraphicalModel<subspace_type,value_type>;

        ConditionedSubmodelBuilder(const gm_type & gm)
        :   m_gm(gm),
            m_sub_num_variables(0),
            m_is_free(m_gm.num_variables(), false),
            m_gm_to_sub_gm(m_gm.num_variables()),
            m_sub_gm_to_gm(m_gm.num_variables()),
            m_factor_fixed_pos(),
            m_factor_fixed_labels()
        {
            const auto max_arity = m_gm.max_arity();

            m_sub_gm_factor_vi.resize(max_arity);
            m_factor_fixed_pos.resize(max_arity);
            m_factor_fixed_labels.resize(max_arity);
        }


        template<class VI_T, class F>
        void condition(std::initializer_list<VI_T> free_vi, const labels_vector_type &  labels, F && f)
        {
            this->condition(free_vi.begin(), free_vi.end(), labels, std::forward<F>(f));
        }


        template<class VI_ITER, class F>
        void condition(VI_ITER free_vi_begin, VI_ITER free_vi_end, const labels_vector_type &  labels, F && f)
        {
            auto subspace = m_gm.space().subspace(free_vi_begin, free_vi_end);
            //std::cout<<"subspace "<<subspace.size()<<" at 0 "<<subspace[0]<<"\n";
            sub_gm_type sub_gm(subspace);
            m_sub_num_variables = sub_gm.num_variables();

            auto svi = 0;
            std::for_each(free_vi_begin, free_vi_end, [&](auto vi){
                m_is_free[vi]=true;
                m_gm_to_sub_gm[vi] = svi;
                m_sub_gm_to_gm[svi] = vi;
                ++svi;
            });

            for(auto && factor : m_gm)
            {
                // * build variable indices of sub factor
                // * get fixed positions
                // * get labels at fixed positions
                auto && vars = factor.variables();
                auto sub_arity = 0;
                auto n_fixed = 0;
                for(auto ai=0; ai<factor.arity(); ++ai)
                {
                    const auto vi = vars[ai];
                    if(m_is_free[vi])
                    {
                        auto sub_vi = m_gm_to_sub_gm[vi];
                        m_sub_gm_factor_vi[sub_arity] = sub_vi;
                        ++sub_arity;
                    }
                    else
                    {
                        m_factor_fixed_pos[n_fixed] = ai;
                        m_factor_fixed_labels[n_fixed] = labels[vi];
                        ++n_fixed;
                    }
                }
                if(sub_arity > 0 && n_fixed > 0)
                {
                    auto subtensor = factor.bind(
                        gsl::span<const std::size_t>(m_factor_fixed_pos.data(), n_fixed),
                        gsl::span<const label_type>(m_factor_fixed_labels.data(), n_fixed)
                    );
                    sub_gm.add_factor(
                        std::move(subtensor),
                        m_sub_gm_factor_vi.begin(),
                        m_sub_gm_factor_vi.begin() + sub_arity
                    );
                }
                else if(sub_arity == factor.arity())
                {
                    auto tensor = factor.tensor();
                    sub_gm.add_factor(
                        tensor,
                        m_sub_gm_factor_vi.begin(),
                        m_sub_gm_factor_vi.begin() + sub_arity
                    );
                }
                // else // exlucded
                // {
                //     std::cout<<"arity "<<factor.arity()<<"\n";
                //     std::cout<<"sub_arity "<<sub_arity<<"\n";
                //     std::cout<<"n_fixed "<<n_fixed<<"\n";
                //     throw std::runtime_error("hups");
                // }

            }
            f(sub_gm);
            // cleanup
            for(auto svi=0; svi<m_sub_num_variables; ++svi)
            {
                m_is_free[m_sub_gm_to_gm[svi]] = true;
            }
        }

    private:

        const gm_type & m_gm;
        std::size_t m_sub_num_variables;
        std::vector<bool> m_is_free;
        std::vector<std::size_t> m_gm_to_sub_gm;
        std::vector<std::size_t> m_sub_gm_to_gm;
        std::vector<std::size_t> m_sub_gm_factor_vi;
        std::vector<std::size_t> m_factor_fixed_pos;
        std::vector<label_type>  m_factor_fixed_labels;
    };

    template<class GM>
    auto conditioned_submodel_builder(const GM & gm)
    {
        return ConditionedSubmodelBuilder<GM>(gm);
    }

}
}