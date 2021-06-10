#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/minimizer/bp.hpp"
#include "opengm/minimizer/utils/label_fuser.hpp"

#include "opengm/factors_of_variables.hpp"

namespace opengm{



template<class GM>
class SelfFusion;


namespace detail_self_fusion{

    template<class GM>
    class ToQuadraticMinimizerVisitor: public MinimizerCallbackBase<GM>{
    public:
        using gm_type = GM;
        using minimizer_base_type =  MinimizerBase<gm_type>;
        using self_fusion_type = SelfFusion<gm_type>;
        using base_type = MinimizerCallbackBase<gm_type>;

        ToQuadraticMinimizerVisitor(self_fusion_type & self_fusion, base_type * outer_visitor)
        :   m_self_fusion(self_fusion),
            m_outer_callback(outer_visitor)
        {
        }

        void begin(minimizer_base_type * minimizer){
            if(m_outer_callback != nullptr){
                return m_outer_callback->begin(&m_self_fusion);
            }
        }
        void end(minimizer_base_type   * minimizer){
            if(m_outer_callback != nullptr){
                return m_outer_callback->end(&m_self_fusion);
            }
        }
        bool operator()(minimizer_base_type  * minimizer){
            if(m_outer_callback != nullptr){
                return m_outer_callback->operator()(&m_self_fusion);
            }
            return true;
        }
    private:

        self_fusion_type & m_self_fusion;
        base_type * m_outer_callback;

    };
}




template<class GM>
class SelfFusion : public MinimizerCrtpBase<GM, SelfFusion<GM> >{
public:

    friend class detail_self_fusion::ToQuadraticMinimizerVisitor<GM>;

    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;
    using base_type::minimize;

    using quadratic_gm_type = typename ToQuadratic<gm_type>::quadratic_gm_type;

    using quadratic_gm_minimizer_factory_base_type = MinimizerFactoryBase<quadratic_gm_type>;
    using quadratic_gm_minimizer_factory_ptr_type = std::shared_ptr<quadratic_gm_minimizer_factory_base_type>;

private:
    using inner_callback_type =  detail_self_fusion::ToQuadraticMinimizerVisitor<gm_type>;

public:






    struct settings_type : public SolverSettingsBase{
        quadratic_gm_minimizer_factory_ptr_type minimizer_factory;
    };


    SelfFusion(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_settings(settings),
        m_energy(),
        m_labels(gm.num_variables(),0),
        m_outer_callback(nullptr),
        m_starting_point_passed(false),
    {
        if(!m_settings.minimizer_factory)
        {
            m_settings.minimizer_factory = std::make_shared<default_minimizer_factory_type>();
        }
    }

    std::string name() const override{
        return "SelfFusion";
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

        // evaluate energy st. callbacks report correct energy
        m_energy = m_gm.evaluate(m_labels);

        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);

        auto quadratic_gm to_quadratic(m_gm);

        auto quadratic_gm_minimizer = m_settings.minimizer_factory->create(quadratic_gm);
        if(m_starting_point_passed && minimizer->can_start_from_starting_point())
        {
            // todo, we need to create the labels for the quadratic gm
            //quadratic_gm_minimizer->set_starting_point(m_labels);
        }

        if(minimizer_callback_base_ptr != nullptr)
        {
            inner_callback_type inner_callback(*this, minimizer_callback_base_ptr);
            quadratic_gm_minimizer->minimize(&inner_callback);
        }
        else
        {
            quadratic_gm_minimizer->minimize();
        }
    }
private:



    const GM & m_gm;
    settings_type m_settings;
    value_type m_energy;
    std::vector<label_type> m_labels;
    bool m_starting_point_passed;
};

}