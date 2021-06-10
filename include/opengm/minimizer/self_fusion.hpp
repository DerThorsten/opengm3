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
    class SelfFusionVisitor: public MinimizerCallbackBase<GM>{
    public:
        using gm_type = GM;
        using minimizer_base_type =  MinimizerBase<gm_type>;
        using self_fusion_type = SelfFusion<gm_type>;
        using base_type = MinimizerCallbackBase<gm_type>;

        SelfFusionVisitor(self_fusion_type & self_fusion, base_type * outer_visitor)
        :   m_self_fusion(self_fusion),
            m_outer_callback(outer_visitor)
        {
        }

        void begin(minimizer_base_type * minimizer){
            // we do nothing here
        }
        void end(minimizer_base_type   * minimizer){
            this->fuse(minimizer);
        }
        bool operator()(minimizer_base_type  * minimizer){
            return this->fuse(minimizer);
        }
    private:

        bool fuse(minimizer_base_type  * minimizer){
            // do the magic
            m_self_fusion.fuse_with(minimizer);

            // call the outer callback
            // (ie the callback pass SelfFusion::minimize)
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
class SelfFusion : public MinimizerCrtpBase<GM, SelfFusion<GM> >{
public:

    friend class detail_self_fusion::SelfFusionVisitor<GM>;

    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;
    using base_type::minimize;

    using label_fuser_type = LabelFuser<gm_type>;
    using label_fuser_settings_type = typename label_fuser_type::settings_type;
    using sub_gm_type = typename label_fuser_type::fuse_gm_type;

private:
    using self_fusion_visitor_type =  detail_self_fusion::SelfFusionVisitor<gm_type>;

public:




    using default_minimizer_factory_type = MinimizerFactory<BeliefPropergation<gm_type>>;
    using minimizer_factory_ptr_type = std::shared_ptr<MinimizerFactoryBase<gm_type>>;
    using fuse_minimizer_factory_ptr_type = typename label_fuser_type::fuse_minimizer_factory_ptr_type;




    struct settings_type : public SolverSettingsBase{
        minimizer_factory_ptr_type minimizer_factory;
        fuse_minimizer_factory_ptr_type fuse_minimizer_factory;
        bool skip_first{true};
    };


    SelfFusion(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_settings(settings),
        m_label_fuser(gm, label_fuser_settings_type{settings.fuse_minimizer_factory}),
        m_energy(),
        m_labels(gm.num_variables(),0),
        m_outer_callback(nullptr),
        m_starting_point_passed(false),
        m_first(true)
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
        auto minimizer = m_settings.minimizer_factory->create(m_gm);
        if(m_starting_point_passed && minimizer->can_start_from_starting_point())
        {
            minimizer->set_starting_point(m_labels);
        }

        self_fusion_visitor_type self_fusion_visitor(*this, minimizer_callback_base_ptr);
        minimizer->minimize(&self_fusion_visitor);

    }
private:

    void fuse_with(base_type  * minimizer){
        const auto  & labels = minimizer->current_labels();
        if(m_first && m_settings.skip_first && !m_starting_point_passed)
        {
            m_labels = labels;
            m_energy = minimizer->current_energy();
        }
        else
        {
            m_energy = m_label_fuser.fuse(labels, m_labels);
        }
        m_first = false;
    };


    const GM & m_gm;
    settings_type m_settings;
    label_fuser_type m_label_fuser;
    value_type m_energy;
    std::vector<label_type> m_labels;

    bool m_starting_point_passed;
    bool m_first;
};

}