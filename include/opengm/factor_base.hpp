#pragma once

#include <array>
#include <gsl-lite/gsl-lite.hpp>

#include "opengm/crtp_base.hpp"
#include "opengm/meta.hpp"


namespace opengm {


    template<class derived>
    class FactorTraits;

    template<class derived>
    class FactorBase : public CrtpBase<derived>{
    public:
        using derived_t = derived;
        using factor_traits_t = FactorTraits<derived_t>;
        using value_type = typename factor_traits_t::value_type;
        using label_type = typename factor_traits_t::label_type;




        value_type operator[](label_type * labels)const{
            return this->derived_cast().operator[](static_cast<const label_type *>(labels));
        }

        std::size_t index(std::size_t vi)const{
            auto && variables = this->derived_cast().variables();
            std::size_t index=0;
            for(auto v : variables){
                if(v == vi){
                    return index;
                }
                ++index;
            }
            return this->derived_cast().arity();
        }

        template<class GM_VAL, class FACTOR_VAL>
        void from_gm(const GM_VAL & gm_val, FACTOR_VAL & factor_val)const{
            auto && variables =  this->derived_cast().variables();
            for(auto i=0; i<variables.size(); ++i){
                factor_val[i] = gm_val[variables[i]];
            }
        }

        template<class ... ARGS, typename  = meta::all_integral<ARGS...>>
        value_type operator()(ARGS && ... args)const{

            std::array<label_type, sizeof ...(ARGS)> labels{
                label_type(std::forward<ARGS>(args)) ...
            };
            return this->derived_cast().operator[](labels.data());
        }

        void copy_corder(value_type * out)const{
            this->derived_cast().tensor()->copy_corder(out);
        }

        void add_values(value_type * out)const{
            this->derived_cast().tensor()->add_values(out);
        }



        void factor_to_variable_messages(
            const value_type ** in_messages,
            value_type ** out_messages)const{
            this->derived_cast().tensor()->factor_to_variable_messages(in_messages, out_messages);
        }


        auto sum_of_shape()const{
            return this->derived_cast().tensor()->sum_of_shape();
        }


        auto bind(
            gsl::span<const std::size_t> positions,
            gsl::span<const label_type> labels
        )const{
            if(this->derived_cast().arity() - positions.size() < 1){
                throw std::runtime_error("cannot bind all variables, at least one variable must be left");
            }
            return this->derived_cast().tensor()->bind(positions, labels);
        }
    };

}