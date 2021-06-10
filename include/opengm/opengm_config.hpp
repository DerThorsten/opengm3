#pragma once
#ifndef OPENGM_OPENGM_CONFIG_HPP
#define OPENGM_OPENGM_CONFIG_HPP



#include "opengm/opengm_version_major.hpp"
#include "opengm/opengm_version_minor.hpp"
#include "opengm/opengm_version_patch.hpp"


namespace opengm {

    using label_type = std::size_t;
    constexpr std::size_t small_vector_arity_size = 5;

} // end namespace opengm


#endif // OPENGM_OPENGM_CONFIG_HPP