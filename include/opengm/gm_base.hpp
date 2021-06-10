#pragma once

#include <array>
#include <vector>
#include "opengm/crtp_base.hpp"

namespace opengm {


    template<class derived>
    class GmTraits;

    template<class derived>
    class GmBase : public CrtpBase<derived>{
    public:
        using derived_t = derived;
        using gm_traits_t = GmTraits<derived_t>;
        using label_type = typename gm_traits_t::label_type;
        using value_type = typename gm_traits_t::value_type;

        using labels_vector_type = std::vector<label_type>;

        auto num_factors()const{
            return std::distance(this->derived_cast().cbegin(), this->derived_cast().cend());
        }

        template<class F>
        void for_each_factor(F && f)const{
            std::size_t fi=0;
            for(auto && factor : this->derived_cast()){
                f(fi, factor);
                ++fi;
            }
        }

        std::size_t size()const{
            return this->derived_cast().space().size();
        }
        std::size_t num_variables()const{
            return this->derived_cast().space().size();
        }
        label_type num_labels(std::size_t vi)const{
            return this->derived_cast().space()[vi];
        }

        std::size_t max_arity() const{
            std::size_t mx = 0;
            for(auto && factor : this->derived_cast()){
                mx = std::max(mx, factor.arity());
            }
            return mx;
        }

        std::size_t arity_upper_bound()const{
            return this->derived_cast().max_arity();
        }




        template<class ITER>
        value_type evaluate_at(ITER labels_begin, ITER labels_end)const{

            if(std::distance(labels_begin, labels_end) != this->derived_cast().num_variables())
            {
                throw std::runtime_error("labels.size != gm.num_variables()");
            }
            for(auto vi=0; vi<this->derived_cast().num_variables(); ++vi)
            {
                if(labels_begin[vi] >= this->derived_cast().num_labels(vi))
                {
                    throw std::runtime_error("label out of bounds");
                }
            }
        }

        template<class CONTAINER>
        value_type evaluate_at(const CONTAINER & labels)const{
            return this->derived_cast().evaluate_at(labels.begin(), labels.end());
        }


        template<class CONTAINER>
        value_type evaluate(const CONTAINER & labels)const{
            return this->derived_cast().evaluate(labels.begin(), labels.end());
        }

        template<class T>
        value_type evaluate(std::initializer_list<T> labels)const{
            return this->derived_cast().evaluate(labels.begin(), labels.end());
        }

        template<class ITER>
        value_type evaluate(ITER labels_begin, ITER labels_end)const{


            std::vector<label_type> label_buffer(this->derived_cast().arity_upper_bound());
            auto energy = value_type(0);
            for(auto && factor : this->derived_cast()){
                auto && variables = factor.variables();
                for( auto i=0; i<variables.size(); ++i){
                    label_buffer[i] = labels_begin[variables[i]];
                }
                energy += factor[label_buffer.data()];
            }
            return energy;

        }
    };



}