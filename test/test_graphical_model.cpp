#include <doctest.h>

#include "opengm/opengm.hpp"
#include "opengm/opengm_config.hpp"


#include "opengm/opengm.hpp"
#include "opengm/space.hpp"
#include "opengm/opengm_config.hpp"
#include "opengm/graphical_model.hpp"



TEST_SUITE_BEGIN("gm");

TEST_CASE("vgm"){

    using value_type = float;
    using value_ilist = std::initializer_list<value_type>;
    using label_type = std::size_t;
    constexpr label_type num_var = 10;
    constexpr label_type num_labels =  2;
    using space_type = opengm::StaticSpace<label_type, num_var, num_labels>;
    using GmType = opengm::GraphicalModel<space_type, value_type>;
    GmType gm;

    auto tid0 = gm.add_function(
        std::make_unique<opengm::Potts2Tensor<value_type>>(2, 1.0)
    );
    CHECK_EQ(tid0, 0);
    auto fid0 = gm.add_factor(tid0, {0,1});
    auto fid1 = gm.add_factor(tid0, {0,1});
    CHECK_EQ(fid0, 0);
    CHECK_EQ(fid1, 1);


    auto tid1 = gm.add_function(
        std::make_unique<opengm::Potts2Tensor<value_type>>(2, 1.0)
    );
    CHECK_EQ(tid1, 1);
    auto fid2 = gm.add_factor(tid1, {0,1});
    auto fid3 = gm.add_factor(tid1, {0,1});
    CHECK_EQ(fid2, 2);
    CHECK_EQ(fid3, 3);

    for(auto && factor : gm){
        factor(1,1);
    }

    gm.evaluate({0,0});



    auto tid2 = gm.add_function(
        std::make_unique<opengm::UnaryTensor<value_type>>(value_ilist{2.0f, 1.0f})
    );
    CHECK_EQ(tid2, 2);
    auto fid4 = gm.add_factor(tid2, {0});
    auto fid5 = gm.add_factor(tid2, {0});


}



TEST_SUITE_END(); // end of testsuite gm
