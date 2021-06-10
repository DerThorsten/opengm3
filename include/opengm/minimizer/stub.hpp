#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/minimizer_base.hpp"

namespace opengm{


template<class GM>
class Stub : public MinimizerCrtpBase<GM, Stub<GM> >{
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


    Stub(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_labels(gm.num_variables(),0),
        m_energy(m_gm.evaluate(m_labels)),
    {
    }

    std::string name() const override{
        return "Stub";
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

    }
private:


    const GM & m_gm;
    std::vector<label_type> m_labels;
    value_type m_energy;

};

template<class GM>
using StubFactory = MinimizerFactory<Stub<GM>>;
}