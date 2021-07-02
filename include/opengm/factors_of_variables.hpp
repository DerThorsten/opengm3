#pragma once


#include <vector>
#include <queue>
#include <set>

namespace opengm{



    template<class GM>
    class FactorsOfVariables : public std::vector<std::vector<std::size_t>>{
    public:
        using base_type =  std::vector<std::vector<std::size_t>>;
        FactorsOfVariables(const GM & gm)
        :   base_type(gm.num_variables()),
            m_gm(gm)
        {
            auto factor_index = std::size_t(0);
            for(auto && factor : gm)
            {
                for(auto var : factor.variables())
                {
                    (*this)[var].push_back(factor_index);
                }
                ++factor_index;
            }
        }

    private:
        const gm_type & m_gm;
    };


    template<class factors_of_variables_type, class F>
    void for_each_neighbour(const factors_of_variables_type & factors_of_variables, const std::size_t var, F && functor)const
    {
        for(auto fi : this->operator[](var))
        {
            for(other_var : m_gm[fi].variables())
            {
                if(var != other_var)
                {
                    functor(var);
                }
            }
        }
    }




    template<class gm_type,class iter_type class factors_of_variables_type>
    auto depth_limited_bfs_search(
        const gm_type & gm, const factors_of_variables_type & factors_of_variables,
        iter_type starting_var_begin,
        iter_type starting_var_end,
        const std::size_t depth_limit
    ){
        while(starting_var_begin != starting_var_end)
        {
            std::set<std::size_t> included;
            included.insert(*starting_var_begin);
            ++starting_var_begin;
        }

        std::queue<std::pair<std::size_t, std::size_t>> queue;
        queue.emplace(starting_var, 0);

        while(!queue.empty())
        {
            if(included.size() >= size_limit)
            {
                break;
            }
            // pop front
            auto [var, d] = queue.front();
            queue.pop();

            if(d < depth_limit)
            {
                for_each_neighbour(factors_of_variables, var, [&](auto & other_var){
                    if(included.find(other_var) == included.end())
                    {
                        queue.emplace(other_var);
                        included.emplace(other_var, d+1);
                    }
                });
            }
        }
        return included;
    }



    namespace detail{
        struct HigherOrderAndUnaryFactorsOfVariablesValueType{
            auto unaries()const{return m_unaries;}
            auto higher_order()const{return m_higher_order;}
            std::vector<std::size_t> m_unaries;
            std::vector<std::size_t> m_higher_order;
        };
    };

    template<class GM>
    class HigherOrderAndUnaryFactorsOfVariables : public std::vector<detail::HigherOrderAndUnaryFactorsOfVariablesValueType>{
    public:
        using base_type = std::vector<detail::HigherOrderAndUnaryFactorsOfVariablesValueType>;
        HigherOrderAndUnaryFactorsOfVariables(const GM & gm)
        : base_type(gm.num_variables()){
            auto factor_index = std::size_t(0);
            for(auto && factor : gm)
            {
                for(auto var : factor.variables())
                {
                    if(factor.arity() == 1){
                        (*this)[var].m_unaries.push_back(factor_index);
                    }
                    else{
                        (*this)[var].m_higher_order.push_back(factor_index);
                    }
                }
                ++factor_index;
            }
        }
    };

}