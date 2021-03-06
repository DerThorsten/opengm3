#pragma once


#include "opengm/factors_of_variables.hpp"

#include <set>

namespace opengm::detail{


    inline int ipow(int base, int exponent)
    {
        int result = 1;
        for (;;)
        {
            if (exponent & 1)
                result *= base;
            exponent >>= 1;
            if (!exponent)
                break;
            base *= base;
        }

        return result;
    }


    // evaluate a labeling for a subset of factors
    template<class GM,class LABELS, class ITER>
    auto evaluate_factors(
        const GM & gm,
        const LABELS & labels,
        ITER factor_indices_begin,
        ITER factor_indices_end,
        std::vector<typename GM::label_type> & buffer
    ){
        auto e = typename GM::value_type(0.0);
        std::for_each(factor_indices_begin, factor_indices_end,  [&](auto fi){
            auto && factor = gm[fi];
            factor.from_gm(labels, buffer);
            e += factor[buffer.data()];
        });
        return e;
    }



    template<class SHAPE, class STATE, class F>
    void for_each_state(const std::size_t size, SHAPE && shape, STATE & state, F && f)
    {
        // initialize labels with zero
        std::fill(state.begin(), state.begin()+size, 0);

        for (;;)
        {

            // do smth with state
            f(state);

            // increment buffered state
            for (size_t vi = size-1; vi >= 0; --vi)
            {
                if (state[vi] <shape[vi] - 1)
                {
                    ++state[vi];
                    break;
                }
                else
                {
                    if(vi != 0)
                    {
                        state[vi] = 0;
                    }
                    else
                    {
                        goto overflow;
                    }
                }
            }

        }
        overflow:
        ;
    }



    template<std::size_t UNIFORM_SHAPE, class STATE, class F>
    void for_each_state(const std::size_t size,  STATE & state, F && f)
    {
        // initialize labels with zero
        std::fill(state.begin(), state.begin()+size, 0);

        for (;;)
        {

            // do smth with state
            f(state);

            // increment buffered state
            for (size_t vi = size-1; vi >= 0; --vi)
            {
                if (state[vi] <UNIFORM_SHAPE- 1)
                {
                    ++state[vi];
                    break;
                }
                else
                {
                    if(vi != 0)
                    {
                        state[vi] = 0;
                    }
                    else
                    {
                        goto overflow;
                    }
                }
            }

        }
        overflow:
        ;
    }
}