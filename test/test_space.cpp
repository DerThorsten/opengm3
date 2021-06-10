#include <doctest.h>

#include "utils.hpp"
#include "opengm/space.hpp"



TEST_SUITE_BEGIN("gm");

TEST_CASE("for_each_state_subset"){

    opengm::ExplicitSpace<std::size_t> space{3, 3, 2, 2};
    using vec_t = std::vector<std::size_t>;
    vec_t labels(space.size());
    labels[0] = 1;
    labels[1] = 1;
    labels[2] = 1;
    labels[3] = 1;
    std::vector<std::size_t> var({1,2});

    auto i=0;
    space.for_each_state(var.begin(), var.end(), labels, [&](auto && lables){
        switch(i){
            case 0:
                CHECK_EQ(labels, vec_t{1,0,0,1});
                break;
            case 1:
                CHECK_EQ(labels, vec_t{1,0,1,1});
                break;
            case 2:
                CHECK_EQ(labels, vec_t{1,1,0,1});
                break;
            case 3:
                CHECK_EQ(labels, vec_t{1,1,1,1});
                break;
            case 4:
                CHECK_EQ(labels, vec_t{1,2,0,1});
                break;
            case 5:
                CHECK_EQ(labels, vec_t{1,2,1,1});
                break;
        }
        ++i;
    });
    CHECK_EQ(i,6);
}


TEST_CASE("for_each_state"){

    opengm::ExplicitSpace<std::size_t> space{2, 3, 2};
    using vec_t = std::vector<std::size_t>;
    vec_t labels(space.size());
    labels[0] = 1;
    labels[1] = 1;
    labels[2] = 1;
    auto i=0;
    space.for_each_state(labels, [&](auto && lables){
        switch(i){
            case 0:
                CHECK_EQ(labels, vec_t{0,0,0}); break;
            case 1:
                CHECK_EQ(labels, vec_t{0,0,1}); break;
            case 2:
                CHECK_EQ(labels, vec_t{0,1,0}); break;
            case 3:
                CHECK_EQ(labels, vec_t{0,1,1}); break;
            case 4:
                CHECK_EQ(labels, vec_t{0,2,0}); break;
            case 5:
                CHECK_EQ(labels, vec_t{0,2,1}); break;
            case 6:
                CHECK_EQ(labels, vec_t{1,0,0}); break;
            case 7:
                CHECK_EQ(labels, vec_t{1,0,1}); break;
            case 8:
                CHECK_EQ(labels, vec_t{1,1,0}); break;
            case 9:
                CHECK_EQ(labels, vec_t{1,1,1}); break;
            case 10:
                CHECK_EQ(labels, vec_t{1,2,0}); break;
            case 11:
                CHECK_EQ(labels, vec_t{1,2,1}); break;
        }
        ++i;
    });
    CHECK_EQ(i,12);


}



TEST_SUITE_END(); // end of testsuite gm
