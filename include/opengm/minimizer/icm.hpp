#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/factors_of_variables.hpp"

namespace opengm{


template<class GM>
class Icm : public MinimizerCrtpBase<GM, Icm<GM> >{
public:
    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;

    using base_type::minimize;

    struct settings_type : public SolverSettingsBase{
        
    };


    Icm(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_labels(gm.num_variables(),0),
        m_factors_of_variables(gm),
        m_value_buffer(m_gm.space().max_num_labels(),0.0),
        m_factor_labels(m_gm.arity_upper_bound(),0),
        m_dirty(),
        m_in_queue(gm.num_variables())
    {
    }

    std::string name() const override{
        return "Icm";
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

        m_energy = m_gm.evaluate(m_labels);

        // start with all variables on the queue
        for(std::size_t vi=0; vi<m_gm.size(); ++vi){
            m_dirty.push(vi);
            m_in_queue[vi] = true;
        }

        while(!m_dirty.empty())
        {
            // get next var
            auto vi = m_dirty.front();
            m_in_queue[vi] = false;
            m_dirty.pop();

            // move the variable vi optimal  where all other variables
            // are conditioned to the then labels as given by "labels"
            const auto changes = this->move_optimal(vi);
            if(changes)
            {
                callback();
                this->set_neighbours_dirty(vi);
            }
        }
    }
private:

    void set_neighbours_dirty(std::size_t vi){
        for(auto && fi : m_factors_of_variables[vi]){
            for(auto other_vi : m_gm[fi].variables())
            {
                if(other_vi != vi && !m_in_queue[other_vi])
                {
                    m_in_queue[other_vi] = true;
                    m_dirty.push(other_vi);
                }
            }
        }
    }

    bool move_optimal(std::size_t vi){

        const auto num_labels = m_gm.num_labels(vi);
        std::fill(m_value_buffer.begin(), m_value_buffer.end(), value_type(0));

        for(auto && vi : m_factors_of_variables[vi]){
            auto && factor = m_gm[vi];
            factor.from_gm(m_labels, m_factor_labels);
            auto vi_pos = factor.index(vi);
            for(auto l=label_type(0); l<num_labels; ++l)
            {
                m_factor_labels[vi_pos] = l;
                m_value_buffer[l] += factor[m_factor_labels.data()];
            }
        }
        const auto old_energy = m_value_buffer[ m_labels[vi]];
        const auto begin_iter = m_value_buffer.begin();
        auto min_label  = std::distance(begin_iter , std::min_element(begin_iter, begin_iter +  num_labels));
        if(min_label != m_labels[vi]){
            m_labels[vi] = min_label;
            const auto new_energy = m_value_buffer[min_label];
            m_energy -= (old_energy - new_energy);
            return true;
        }
        return false;
    };



    const GM & m_gm;
    value_type m_energy;
    std::vector<label_type> m_labels;
    FactorsOfVariables<GM> m_factors_of_variables;

    std::vector<value_type> m_value_buffer;
    std::vector<label_type> m_factor_labels;

    std::queue<std::size_t> m_dirty;
    std::vector<bool> m_in_queue;
};

template<class GM>
using IcmFactory = MinimizerFactory<Icm<GM>>;
}