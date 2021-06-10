#pragma once


#include "opengm/factors_of_variables.hpp"
#include "opengm/utils.hpp"

#include <set>
#include <iostream>

namespace opengm{






template<class GM>
class Movemaker{
public:
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename GM::labels_vector_type;

    Movemaker(const GM & gm)
    :   m_gm(gm),
        m_factors_of_varaibles(gm),
        m_labels(gm.num_variables(), 0),
        m_labels_buffer(gm.num_variables(), 0),
        m_labels_buffer2(gm.num_variables(), 0),
        m_factor_labels(gm.max_arity(),0)
    {
        m_energy =  m_gm.evaluate(m_labels.begin(), m_labels.end());
    }


    template<class ITER, class F>
    bool move_brute_force_optimal(ITER vars_begin, ITER vars_end, F && f){

        //
        auto factors = this->factors_of_variables(vars_begin, vars_end);
        auto current_e = this->evaluate_factors(m_labels, factors.begin(), factors.end());
        auto best_e = current_e;

        // labels before
        std::for_each(vars_begin, vars_end, [&](auto var){
            m_labels_buffer2[var] = m_labels[var];
        });

        // iterate over all states of the subset of variables
        m_gm.space().for_each_state(vars_begin, vars_end, m_labels_buffer, [&](auto && labels){

            // get energy of this labeling / configuration
            auto e = this->evaluate_factors(labels, factors.begin(), factors.end());

            // found new best energy?
            if(e < best_e)
            {
                std::for_each(vars_begin, vars_end, [&](auto var){
                    m_labels[var] = labels[var];
                });
                best_e = e;
            }

        });

        bool any_change = false;
        // make sure buffers are in sync
        std::for_each(vars_begin, vars_end, [&](auto var){
            if(m_labels[var] != m_labels_buffer2[var])
            {
                any_change = true;
            }
            m_labels_buffer[var] = m_labels[var];
        });

        if(any_change){
            f(factors);
            m_energy -= (current_e - best_e);
            return true;
        }
        else{
            return false;
        }
    }

    template<class ITER>
    bool move_brute_force_optimal(ITER vars_begin, ITER vars_end){
        return this->move_brute_force_optimal(vars_begin, vars_end, [&](auto && factors){

        });
    }

    template<class ITER>
    void set_labels(ITER labels_begin, ITER labels_end){
        m_labels.assign(labels_begin, labels_end);
        m_labels_buffer.assign(labels_begin, labels_end);
        m_energy = m_gm.evaluate(m_labels);
    }

    const auto & factors_of_variables()const{
        return m_factors_of_varaibles;
    }
    const auto & labels()const{
        return m_labels;
    }
    auto energy()const{
        return m_energy;
    }
private:

    // get the set of factors for a set of variables
    template<class ITER>
    auto factors_of_variables(ITER vars_begin, ITER vars_end){

        std::set<size_t> factors;
        std::for_each(vars_begin, vars_end, [&](auto var)
        {
            for(auto f : m_factors_of_varaibles[var])
            {
                factors.insert(f);
            }
        });
        return factors;
    }

    // evaluate a labeling for a subset of factors
    template<class ITER>
    auto evaluate_factors(const labels_vector_type & labels, ITER factor_indices_begin, ITER factor_indices_end){
        return detail::evaluate_factors(m_gm, labels, factor_indices_begin, factor_indices_end, m_factor_labels);
    }







    const GM & m_gm;
    value_type m_energy;
    FactorsOfVariables<GM> m_factors_of_varaibles;
    labels_vector_type m_labels;
    labels_vector_type m_labels_buffer;
    labels_vector_type m_labels_buffer2;
    std::vector<label_type> m_factor_labels;

};

}