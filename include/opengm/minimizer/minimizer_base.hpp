#pragma once

#include <queue>
#include <algorithm>
#include <string>
#include <iostream>
#include <memory>

#include "opengm/from_gm_factory.hpp"
#include "opengm/crtp_base.hpp"

namespace opengm{

class SolverSettingsBase{
public:
};



template<class GM>
class MinimizerBase;


template<class GM>
class MinimizerCallbackBase{
public:
    using gm_type = GM;
    using minimizer_base_type =  MinimizerBase<GM>;


    virtual void begin(minimizer_base_type * minimizer_base) = 0;
    virtual void end(minimizer_base_type   * minimizer_base) = 0;
    virtual bool operator()(minimizer_base_type  * minimizer_base) = 0;
private:
};


template<class GM>
class VerboseMinimizerCallback : public MinimizerCallbackBase<GM>{
public:
    using minimizer_base_type =  MinimizerBase<GM>;

    VerboseMinimizerCallback(const std::size_t visit_nth = 1, const bool same_line=true)
    :   m_iteration(0),
        m_visit_nth(visit_nth),
        m_same_line(same_line)
    {
    }

    void begin(minimizer_base_type * minimizer) override{
        std::cout<<"begin: current " <<minimizer->current_energy()<<" best "<<minimizer->best_energy()<<"\n";
    }
    void end(minimizer_base_type   * minimizer) override{
        if(m_same_line){
            std::cout<<"\n";
        }
        std::cout<<"end:  current " <<minimizer->current_energy()<<" best "<<minimizer->best_energy()<<"\n";
    }
    constexpr bool operator()(minimizer_base_type  * minimizer) override{
        if(m_iteration == 0 || (m_iteration % m_visit_nth == 0))
        {
            if(m_same_line){
                std::cout<<"\r"<<m_iteration<<": current " <<minimizer->current_energy()<<" best "<<minimizer->best_energy()<<std::flush;
            }
            else{
                std::cout<<m_iteration<<": current " <<minimizer->current_energy()<<" best "<<minimizer->best_energy()<<"\n";
            }
        }
        ++m_iteration;
        return true;
    }
private:
    std::size_t m_iteration;
    std::size_t m_visit_nth;
    bool m_same_line;
};





template<class GM>
class CallbackWrapper{
public:

    using gm_type = GM;
    using minimizer_base_type =  MinimizerBase<GM>;
    using minimizer_callback_base_type =  MinimizerCallbackBase<gm_type>;
    using minimizer_callback_base_ptr_type =  minimizer_callback_base_type *;
    using minimizer_callback_base_unique_ptr_type =  std::unique_ptr<minimizer_callback_base_type>;

    CallbackWrapper(minimizer_base_type * minimizer, minimizer_callback_base_ptr_type minimizer_callback_base_ptr)
    :   m_minimizer(minimizer),
        m_minimizer_callback_base_ptr(std::move(minimizer_callback_base_ptr)){

        if(m_minimizer_callback_base_ptr != nullptr){
            m_minimizer_callback_base_ptr->begin(m_minimizer);
        }
    }
    ~CallbackWrapper(){
        if(m_minimizer_callback_base_ptr != nullptr){
            m_minimizer_callback_base_ptr->end(m_minimizer);
        }
    }

    bool operator()() {
        if(m_minimizer_callback_base_ptr != nullptr){
            return m_minimizer_callback_base_ptr->operator()(m_minimizer);
        }
        return true;
    }
    operator bool() const {
        return m_minimizer_callback_base_ptr != nullptr
    }

private:
    minimizer_base_type * m_minimizer;
    minimizer_callback_base_ptr_type m_minimizer_callback_base_ptr;
};


template<class GM>
auto callback_wrapper( MinimizerBase<GM> * minimizer, MinimizerCallbackBase<GM> * minimizer_callback_base_ptr){
    return CallbackWrapper(minimizer, minimizer_callback_base_ptr);
}

// template<class GM>
// class MinimizerFactoryBase;


template<class GM>
class MinimizerBase{
public:
    using gm_type = GM;

    //using minimizer_factory_base = MinimizerFactoryBase<gm_type>;
    using minimizer_callback_base_type =  MinimizerCallbackBase<gm_type>;
    using minimizer_callback_base_ptr_type =  minimizer_callback_base_type *;
    using minimizer_callback_base_unique_ptr_type =  std::unique_ptr<minimizer_callback_base_type>;
    using value_type = typename gm_type::value_type;
    using label_type = typename gm_type::label_type;


    using labels_vector_type = std::vector<label_type>;

    virtual void minimize(minimizer_callback_base_unique_ptr_type callback){
        this->minimize(callback.get());
    }
    virtual void minimize(minimizer_callback_base_ptr_type) = 0;
    virtual void minimize(){
        this->minimize(minimizer_callback_base_ptr_type());
    }
    virtual bool can_start_from_starting_point() = 0;
    virtual std::string name() const = 0;

    virtual void set_starting_point(const label_type * labels){
    }
    virtual void set_starting_point(const labels_vector_type & labels){
        this->set_starting_point(labels.data());
    }

    virtual const gm_type & gm() const = 0;
    virtual const labels_vector_type & best_labels() = 0;
    virtual const labels_vector_type & current_labels() = 0;
    virtual value_type best_energy() {
        return this->gm().evaluate(this->best_labels());
    }
    virtual value_type current_energy() {
        return this->gm().evaluate(this->current_labels());
    }
};



template<class GM, class DERIVED>
class MinimizerCrtpBase : public MinimizerBase<GM>, public CrtpBase<DERIVED>
{
public:
};




template<class gm_type>
using MinimizerFactoryBase =  detail::FromGmFactoryBase<
    MinimizerBase<gm_type>,
    gm_type
>;


template<class MINIMIZER>
using MinimizerFactory = detail::FromGmFactory<
    MINIMIZER,
    MinimizerBase<typename MINIMIZER::gm_type>,
    typename MINIMIZER::gm_type,
    typename MINIMIZER::settings_type
>;


template<class MINIMIZER, class F>
auto make_shared_factory(F && f){
    using settings_type = typename MINIMIZER::settings_type;
    settings_type settings;
    f(settings);
    return std::make_shared<MinimizerFactory<MINIMIZER>>(settings);
}

template<class MINIMIZER>
auto make_shared_factory(){
    return std::make_shared<MinimizerFactory<MINIMIZER>>();
}


}