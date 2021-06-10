#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <utility>


#include "opengm/meta.hpp"
#include "opengm/factor_base.hpp"
#include "opengm/factors.hpp"
#include "opengm/space.hpp"
#include "opengm/gm_base.hpp"
#include "opengm/tensors.hpp"

namespace opengm{


  

    template<class SPACE, class T, class DERIVED >
    class TensorViewGm;

    template<class SPACE, class T, class DERIVED>
    class GmTraits<TensorViewGm<SPACE,T, DERIVED>>{
    public:
        using space_type = SPACE;
        using label_type = typename SpaceTraits<space_type>::label_type;
        using value_type = T;
    };




    template<class SPACE, class T, class DERIVED =  meta::null_type >
    class TensorViewGm : public GmBase<meta::if_not_null_type_t<DERIVED, TensorViewGm<SPACE, T, DERIVED>>>
    {
    public:

        using space_type = SPACE;
        using tensor_type = TensorBase<T>;

        template<class ... ARGS>
        TensorViewGm(ARGS && ... args)
        :   m_space(std::forward<ARGS>(args)...),
            m_factors()
        {

        }

        TensorViewGm(space_type && space)
        :   m_space(space),
            m_factors()
        {

        }

        template<class ITER>
        auto add_factor(tensor_type * tensor, ITER var_begin, ITER var_end){
            const auto fid = m_factors.size();
            m_factors.emplace_back(tensor, var_begin, var_end);
            return fid;
        }


        template<class ITER>
        auto add_factor(const tensor_type * tensor, ITER var_begin, ITER var_end){
            const auto fid = m_factors.size();
            m_factors.emplace_back(tensor, var_begin, var_end);
            return fid;
        }

        auto cbegin() const{
            return m_factors.cbegin();
        }
        auto cend() const{
            return m_factors.cend();
        }

        auto begin() const{
            return m_factors.begin();
        }
        auto end() const{
            return m_factors.end();
        }


        const auto & space()const{
            return m_space;
        }
        auto & space(){
            return m_space;
        }

        const auto & operator[](const std::size_t fi)const{
            return m_factors[fi];
        }
        void clear(){
            m_factors.clear();
        }
    private:

        space_type m_space;
        std::vector<VFactor<T>> m_factors;
    };
}