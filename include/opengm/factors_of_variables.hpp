#pragma once


#include <vector>

namespace opengm{



    template<class GM>
    class FactorsOfVariables : public std::vector<std::vector<std::size_t>>{
    public:
        using base_type =  std::vector<std::vector<std::size_t>>;
        FactorsOfVariables(const GM & gm)
        : base_type(gm.num_variables()){
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
    };



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