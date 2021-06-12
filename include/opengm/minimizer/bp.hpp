#pragma once

#include <queue>
#include <algorithm>
#include <map>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/arity_vector.hpp"

namespace opengm{





namespace detail{

    template<class T>
    struct Msg{
        T * msg_;
        T * oMsg_;
    };



    template<class GM>
    class MessageStoring{
    public:
        typedef GM gm_type;
        using value_type = typename gm_type::value_type;
        using factors_of_variables_type =  HigherOrderAndUnaryFactorsOfVariables<gm_type>;
        MessageStoring(const gm_type & gm, const factors_of_variables_type & factors_of_variables)
        :   m_gm(gm),
            m_factors_of_variables(factors_of_variables),
            m_msg_ptrs(),
            m_fac_to_var_offset(gm.num_factors()),
            m_var_to_fac_offset(gm.num_variables())
        {

            uint64_t num_fac_to_var_msg = 0;
            uint64_t msg_value_space = 0;
            // count messages
            for(auto && factor : m_gm){
                const auto arity = factor.arity();
                if(arity>1){
                    num_fac_to_var_msg +=arity;
                    msg_value_space += factor.sum_of_shape();
                }
            }

            m_msg_storage.resize(2*msg_value_space,0.0);
            m_msg_ptrs.reserve(2*num_fac_to_var_msg);

            // setup offsets
            uint64_t var_offset = 0;
            uint64_t msg_index = 0;

            for(auto vi=0; vi<m_gm.num_variables(); ++vi){
                const auto nl = m_gm.num_labels(vi);
                const auto & factors = m_factors_of_variables[vi];
                m_var_to_fac_offset[vi] = msg_index;
                for(auto fi : factors.higher_order()){
                    const auto factor = m_gm[fi];

                    m_msg_ptrs.push_back(Msg<value_type>{m_msg_storage.data()+ var_offset , NULL});
                    var_offset += nl;
                    ++msg_index;
                }
            }

            m_gm.for_each_factor([&](auto fi, auto && factor){
                const auto arity = factor.arity();
                if(arity > 1){

                    auto && variables = factor.variables();
                    // set var_offset
                    for(auto a=0; a<arity; ++a){


                        const auto vi = variables[a];
                        const auto nl = m_gm.num_labels(vi);
                        auto && hfacs = m_factors_of_variables[vi].higher_order();

                        // find out at which position fDesc is in hfacs
                        auto fi_iter = std::find(hfacs.begin(), hfacs.end(), fi);
                        const auto pos = std::distance(hfacs.begin(), fi_iter);

                        auto  & varToFacMsgHolder = m_msg_ptrs[m_var_to_fac_offset[vi]+pos];

                        const auto facToVarMsgHolder  = Msg<value_type>{m_msg_storage.data()+ var_offset, varToFacMsgHolder.msg_};
                        varToFacMsgHolder.oMsg_ = facToVarMsgHolder.msg_;
                        m_msg_ptrs.push_back(facToVarMsgHolder);
                        var_offset += nl;
                    }
                    m_fac_to_var_offset[fi] = msg_index;
                    msg_index += arity;
                }
            });
        }

        value_type * facToVarMsg(const std::size_t fi, const size_t mi){
            const uint64_t offset = m_fac_to_var_offset[fi];
            return m_msg_ptrs[offset + mi].msg_;
        }

        value_type * oppToFacToVarMsg(const std::size_t fi, const size_t mi){
            const uint64_t offset = m_fac_to_var_offset[fi];
            return m_msg_ptrs[offset + mi].oMsg_;
        }

        value_type * varToFacMsg(const std::size_t vi, const size_t mi){
            const uint64_t offset = m_var_to_fac_offset[vi];
            return m_msg_ptrs[offset + mi].msg_;
        }

        value_type * oppToVarToFacMsg(const std::size_t vi, const size_t mi){
            const uint64_t offset = m_var_to_fac_offset[vi];
            return m_msg_ptrs[offset + mi].oMsg_;
        }


        uint64_t nMsg()const{
            return m_msg_ptrs.size();
        }

    private:
        const gm_type & m_gm;
        const factors_of_variables_type & m_factors_of_variables;
        std::vector<Msg<value_type>> m_msg_ptrs;

        std::vector<std::size_t> m_fac_to_var_offset;
        std::vector<std::size_t> m_var_to_fac_offset;
        std::vector<value_type> m_msg_storage;
        // uint64_t lastVarOffset_;
    };











}





template<class GM>
class BeliefPropergation : public MinimizerCrtpBase<GM, BeliefPropergation<GM> >{
    //public MinimizerBase<GM>{
public:


    using gm_type = GM;
    using base_type = MinimizerBase<gm_type>;
    using factors_of_variables_type =  HigherOrderAndUnaryFactorsOfVariables<gm_type>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;

    using base_type::minimize;

    struct Settings : public SolverSettingsBase{
        std::size_t num_iterations{10000};
        value_type damping{0.9};
        value_type convergence{5e-7};
    };

    using settings_type = Settings;


    struct MsgBase{
        MsgBase(std::size_t v, std::size_t f, std::size_t nl)
        :   vi(v),
            fi(v),
            values(nl, value_type(0)){
        }
        auto data() {
            return values.data();
        }

        std::size_t vi;
        std::size_t fi;
        std::vector<value_type> values;
    };

    struct VarToFacMsg : public MsgBase{
        using MsgBase::MsgBase;
    };
    struct FacToVarMsg : public MsgBase{
        using MsgBase::MsgBase;
    };
 

    BeliefPropergation(const GM & gm, const Settings & settings = Settings())
    :   m_gm(gm),
        m_settings(settings),
        m_factors_of_variables(gm),
        m_msg(gm, m_factors_of_variables),
        m_current_energy(),
        m_best_energy(),
        m_current_labels(gm.num_variables(), 0),
        m_best_labels(gm.num_variables(),0),
        sMsgBuffer_(gm.space().max_num_labels())
    {
        m_current_energy =  m_gm.evaluate(m_current_labels);
        m_best_energy = m_current_energy;

    }

    std::string name() const override{
        return "BeliefPropergation";
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


        for(auto iteration=0; iteration<m_settings.num_iterations; ++iteration)
        {
            this->sendAllFacToVar();
            auto eps = this->sendAllVarToFac();
            m_current_energy =  m_gm.evaluate(m_current_labels);
            callback();
            if(m_current_energy < m_best_energy)
            {
                m_best_energy = m_current_energy;
                m_best_labels = m_current_labels;
            }

            if( eps < m_settings.convergence)
            {
                break;
            }
        };

    }

    auto sendAllVarToFac(){
        auto eps = value_type(0);
        for(auto vi=0; vi<m_gm.num_variables(); ++vi)
        {
            eps +=  this->sendVarToFac(vi);
        }
        eps /= m_msg.nMsg();
        return eps;
    }

    void sendAllFacToVar(){
        for(auto fi=0; fi<m_gm.num_factors(); ++fi)
        {
            this->sendFacToVar(fi);
        }
    }


    void sendFacToVar(const std::size_t fi){
        auto && factor = m_gm[fi];
        const auto arity  = factor.arity();

        arity_vector<value_type *>       facToVar(arity);
        arity_vector<const value_type *> varToFac(arity);
        for(auto i=0; i<arity; ++i){
            facToVar[i] = m_msg.facToVarMsg(fi, i);
            varToFac[i] = m_msg.oppToFacToVarMsg(fi, i);
        }
        factor.factor_to_variable_messages(varToFac.data(), facToVar.data());
    }


    value_type sendVarToFac(const std::size_t vi){

        auto msg_squared_diff = value_type(0.0);

        auto buffer = sMsgBuffer_.data();

        // how many labels
        const auto num_labels = m_gm.num_labels(vi);

        // factors for this variable
        auto && unaries = m_factors_of_variables[vi].unaries();
        auto && higher_order = m_factors_of_variables[vi].higher_order();

        if(unaries.size() + higher_order.size() > 0)
        {
            // initialize buffer
            std::fill(buffer, buffer + num_labels, 0.0);

            // add unaries to buffer
            for(auto fi : unaries)
            {
                m_gm[fi].add_values(buffer);
            }

            // higher order factors
            for(auto hoi=0; hoi<higher_order.size(); ++hoi)
            {
                const auto fac_to_var = m_msg.oppToVarToFacMsg(vi, hoi);
                for(label_type l=0; l<num_labels; ++l)
                {
                    buffer[l] +=fac_to_var[l];
                }
            }

            // all msg are summed up now therefore buffer is
            // the actual belief vector
            m_current_labels[vi] = std::distance(buffer, std::min_element(buffer,buffer+num_labels));


            // compute outgoing messages
            // - we subtract a single varToFac msg 
            //   from belief vector 
            for(size_t hoi=0; hoi<higher_order.size(); ++hoi){
                // get the factor msg
                const auto fac_to_var = m_msg.oppToVarToFacMsg(vi,hoi);
                auto var_to_fac = m_msg.varToFacMsg(vi,hoi);

                // calculate the mean of the new message
                value_type mean = 0.0;
                for(label_type l=0; l<num_labels; ++l){
                    mean += (buffer[l] - fac_to_var[l]);
                }
                mean /= num_labels;


              
                // substract the mean and damp
                for(label_type l=0; l<num_labels; ++l){
                    const value_type oldValue = var_to_fac[l];
                    const value_type newUndampedValue = (buffer[l] - fac_to_var[l])- mean;
                    const value_type newDampedValue = m_settings.damping * oldValue + (1.0 - m_settings.damping)*newUndampedValue;
                    var_to_fac[l] = newDampedValue;
                    // convergence accumulation
                    const value_type diff = oldValue-newDampedValue;
                    msg_squared_diff += diff*diff;
                }
            }
        }
        return msg_squared_diff;
    }
private:




    const gm_type & m_gm;
    Settings m_settings;
    factors_of_variables_type m_factors_of_variables;
    detail::MessageStoring<GM> m_msg;
    value_type m_current_energy;
    value_type m_best_energy;
    labels_vector_type m_current_labels;
    labels_vector_type m_best_labels;
    std::vector<value_type> sMsgBuffer_;


};



template<class GM>
using BeliefePropergationFactory = MinimizerFactory<BeliefPropergation<GM>>;




}