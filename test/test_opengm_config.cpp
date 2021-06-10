#include <doctest.h>

#include "opengm/opengm.hpp"
#include "opengm/opengm_config.hpp"



TEST_SUITE_BEGIN("core");

TEST_CASE("check version"){

    #ifndef OPENGM_VERSION_MAJOR
        #error "OPENGM_VERSION_MAJOR is undefined"
    #endif

    #ifndef OPENGM_VERSION_MINOR
        #error "OPENGM_VERSION_MINOR is undefined"
    #endif


    #ifndef OPENGM_VERSION_PATCH
        #error "OPENGM_VERSION_PATCH is undefined"
    #endif

    CHECK_EQ(OPENGM_VERSION_MAJOR , 0);
    CHECK_EQ(OPENGM_VERSION_MINOR , 1);
    CHECK_EQ(OPENGM_VERSION_PATCH , 0);
}



TEST_SUITE_END(); // end of testsuite core
