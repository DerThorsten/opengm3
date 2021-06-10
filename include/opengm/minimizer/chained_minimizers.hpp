#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/minimizer/bp.hpp"
#include "opengm/minimizer/utils/label_fuser.hpp"

#include "opengm/factors_of_variables.hpp"

namespace opengm{



template<class GM>
class ChainedMinimizers;


namespace detail_self_fusion{

    template<class GM>
    class ChaindedMinimizersVisitor: public MinimizerCallbackBase<GM>{
    public:
        using gm_type = GM;
        using minimizer_base_type =  MinimizerBase<gm_type>;
        using self_fusion_type = ChainedMinimizers<gm_type>;
        using base_type = MinimizerCallbackBase<gm_type>;

        ChaindedMinimizersVisitor(self_fusion_type & self_fusion, base_type * outer_visitor)
        :   m_self_fusion(self_fusion),
            m_outer_callback(outer_visitor)
        {
        }

        void begin(minimizer_base_type * minimizer){
            // we do nothing here
        }
        void end(minimizer_base_type   * minimizer){
            this->call_outer_callack(minimizer);
        }
        bool operator()(minimizer_base_type  * minimizer){
            return this->call_outer_callack(minimizer);
        }
    private:

        bool call_outer_callack(minimizer_base_type  * minimizer){
            // call the outer callback
            // (ie the callback pass ChainedMinimizers::minimize)
            if(m_outer_callback != nullptr){
                return m_outer_callback->operator()(&m_self_fusion);
            }
            else{
                return true;
            }

        }
        self_fusion_type & m_self_fusion;
        base_type * m_outer_callback;

    };
}




template<class GM>
class ChainedMinimizers : public MinimizerCrtpBase<GM, ChainedMinimizers<GM> >{
public:

    friend class detail_self_fusion::ChaindedMinimizersVisitor<GM>;

    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;
    using base_type::minimize;


    using chained_minimizer_visitor_type =  detail_self_fusion::ChaindedMinimizersVisitor<gm_type>;

public:


    using minimizer_factory_base_type = MinimizerFactoryBase<gm_type>;
    using minimizer_factory_ptr_type = std::shared_ptr<minimizer_factory_base_type>;


    struct settings_type : public SolverSettingsBase{
        std::vector<minimizer_factory_ptr_type> minimizer_factories;
    };


    ChainedMinimizers(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_settings(settings),
        m_energy(),
        m_labels(gm.num_variables(),0),
        m_outer_callback(nullptr),
        m_starting_point_passed(false)
    {
    }

    std::string name() const override{
        return "ChainedMinimizers";
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
        m_starting_point_passed = true;
        m_labels.assign(labels, labels + m_gm.num_variables());
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override{

        // store callback ptr st we can call it from the chained_minimizer_visitor
        m_outer_callback = minimizer_callback_base_ptr;

        // evaluate energy st. callbacks report correct energy
        m_energy = m_gm.evaluate(m_labels);

        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);

        chained_minimizer_visitor_type chained_minimizer_visitor(*this, minimizer_callback_base_ptr);

        bool first = true;
        for(auto minimizer_factory : m_settings.minimizer_factories)
        {
            auto minimizer = minimizer_factory->create(m_gm);
            if(!first ||  m_starting_point_passed && minimizer->can_start_from_starting_point())
            {
                minimizer->set_starting_point(m_labels);
            }
            minimizer->minimize(&chained_minimizer_visitor);
            const auto minimizer_energy = minimizer->best_energy();
            if(minimizer_energy < m_energy)
            {
                m_labels = minimizer->best_labels();
                m_energy = minimizer_energy;
            }
            first = false;
        }

        m_starting_point_passed = false;
        m_outer_callback = nullptr;
    }
private:


    const GM & m_gm;
    settings_type m_settings;
    value_type m_energy;
    std::vector<label_type> m_labels;

    minimizer_callback_base_ptr_type m_outer_callback;

    bool m_starting_point_passed;
};

}