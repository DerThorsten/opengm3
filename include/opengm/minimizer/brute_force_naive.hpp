#pragma once

#include <queue>
#include <algorithm>
#include <map>

#include "opengm/minimizer/minimizer_base.hpp"

namespace opengm{


template<class GM>
class BruteForceNaive : public MinimizerCrtpBase<GM, BruteForceNaive<GM> >{
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

    BruteForceNaive(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_current_labels(gm.num_variables(), 0)
    {

    }

    std::string name() const override{
        return "BruteForceNaive";
    }
    const gm_type & gm() const override{
        return m_gm;
    }

    const labels_vector_type & best_labels()override{
        return m_best_labels;
    }
    const labels_vector_type & current_labels() override{
        return m_current_labels;
    }
    value_type best_energy()  override{
        return m_best_energy;
    }
    value_type current_energy()  override{
        return m_current_energy;
    }

    bool can_start_from_starting_point() override{
        return false;
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override{
        m_current_energy =  m_gm.evaluate(m_current_labels);
        m_best_labels = m_current_labels;
        m_best_energy = m_current_energy;

        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);

        // naive
        m_gm.space().for_each_state(m_current_labels, [&](auto && lables){
            m_current_energy = m_gm.evaluate(m_current_labels);
            if(m_current_energy < m_best_energy)
            {
                callback();
                m_best_energy = m_current_energy;
                m_best_labels = m_current_labels;
            }
        });
    }
private:




    const GM & m_gm;
    value_type m_current_energy;
    value_type m_best_energy;
    labels_vector_type m_current_labels;
    labels_vector_type m_best_labels;
};


template<class GM>
using BruteForceNaiveFactory = MinimizerFactory<BruteForceNaive<GM>>;

}