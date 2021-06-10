#pragma once

#include <vector>
#include <memory>
#include <utility>


#include "opengm/factor_base.hpp"
#include "opengm/factors.hpp"
#include "opengm/tensors.hpp"

namespace opengm{
    template<class T>
    class VFactor;

    template<class T>
    class FactorTraits<VFactor<T>>{
    public:
        using value_type = T;
        using label_type = std::size_t;
    };


    template<class T>
    class VFactor : public FactorBase<VFactor<T>>{
    public:
        using base_type = FactorBase<VFactor<T>>;
        using base_type::operator();
        using base_type::operator[];
        using label_type = std::size_t;
        using value_type = T;

        template<class ITER>
        VFactor(const TensorBase<T> *  tensor, ITER var_begin, ITER var_end)
        :   m_tensor(tensor),
            m_variables(var_begin, var_end){
        }
        const auto & variables()const{
            return m_variables;
        }

        virtual std::size_t arity()const{
            return m_tensor->arity();
        }
        virtual std::size_t shape(const std::size_t i) const{
            return m_tensor->shape(i);
        }
        value_type operator[](const label_type * labels)const{
            return m_tensor->operator[](labels);
        }

        const TensorBase<T> * tensor() const{
            return m_tensor;
        }


    private:
        const TensorBase<T> * m_tensor;
        std::vector<std::size_t> m_variables;
    };
}