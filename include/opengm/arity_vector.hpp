#pragma once

#include <boost/container/small_vector.hpp>


#include "opengm/opengm_config.hpp"

namespace opengm {



    template<class T>
    using arity_vector = boost::container::small_vector<T, small_vector_arity_size>;

} // end namespace opengm
