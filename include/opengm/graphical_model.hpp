#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <utility>

#include "opengm/factor_base.hpp"
#include "opengm/factors.hpp"
#include "opengm/space.hpp"
#include "opengm/gm_base.hpp"
#include "opengm/tensors.hpp"
#include "opengm/tensor_view_gm.hpp"

namespace opengm{


    

    template<class SPACE, class T>
    class GraphicalModel;

    template<class SPACE, class T>
    class GmTraits<GraphicalModel<SPACE,T>>{
    public:
        using space_type = SPACE;
        using label_type = typename SpaceTraits<space_type>::label_type;
        using value_type = T;
    };



    template<class SPACE, class T>
    class GraphicalModel : public TensorViewGm<SPACE, T, GraphicalModel<SPACE, T>>{
    public:

        using base_type =  TensorViewGm<SPACE, T, GraphicalModel<SPACE, T>>;
        using space_type = SPACE;
        using virtual_tensor_base_type = TensorBase<T>;
        using tensor_type = virtual_tensor_base_type;
        using unique_tensor_ptr = std::unique_ptr<tensor_type>;

        template<class ... ARGS>
        GraphicalModel(ARGS && ... args)
        :   base_type(std::forward<ARGS>(args)...),
            m_tensors()
        {

        }

        GraphicalModel(space_type && space)
        :   base_type(std::forward<space_type>(space)),
            m_tensors()
        {

        }

        auto add_tensor(unique_tensor_ptr tensor){
            const auto tid = m_tensors.size();
            m_tensors.push_back(std::move(tensor));
            return  tid;
        }

        auto add_function(unique_tensor_ptr tensor){
            const auto tid = m_tensors.size();
            m_tensors.push_back(std::move(tensor));
            return  tid;
        }

        template<class ITER>
        auto add_factor(std::size_t tid, ITER var_begin, ITER var_end){
            return base_type::add_factor(m_tensors[tid].get(), var_begin, var_end);
        }
        template<class VAR_T>
        auto add_factor(std::size_t tid, std::initializer_list<VAR_T> vars){
            return this->add_factor(tid, vars.begin(), vars.end());
        }

        template<class VAR_T>
        auto add_factor(unique_tensor_ptr tensor,std::initializer_list<VAR_T> vars){
            return this->add_factor(std::move(tensor), vars.begin(), vars.end());
        }

        auto add_unary_factor(std::size_t tid, std::size_t vi){
            return this->add_factor(tid, &vi, &vi+1);
        }

        auto add_unary_factor(unique_tensor_ptr tensor, std::size_t vi){
            auto function_id = add_tensor(std::move(tensor));
            return this->add_unary_factor(function_id, vi);
        }


        template<class ITER>
        auto add_factor(unique_tensor_ptr tensor, ITER var_begin, ITER var_end){
            auto function_id = add_tensor(std::move(tensor));
            return this->add_factor(function_id, var_begin, var_end);
        }

        template<class ITER>
        auto add_factor(const tensor_type * tensor, ITER var_begin, ITER var_end){
            return base_type::add_factor(tensor, var_begin, var_end);
        }


        void clear(){
            base_type::clear();
            m_tensors.clear();
        }
    private:
        std::vector<unique_tensor_ptr> m_tensors;
    };
}