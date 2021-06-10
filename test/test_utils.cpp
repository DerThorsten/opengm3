#include <doctest.h>
#include "utils.hpp"

#include "opengm/space.hpp"
#include "opengm/utils.hpp"
#include "opengm/toy_models.hpp"



TEST_SUITE_BEGIN("gm");

TEST_CASE("utils"){

    auto n_variables = 4;
    auto n_labels =  4;
    auto seed = 42;
    auto gm = opengm::RandomPottsChain(n_variables, n_labels, seed)();
    using gm_type = std::decay_t<decltype(gm)>;
    using label_type = typename gm_type::label_type;
    using vec_t = std::vector<std::size_t>;
    CHECK_EQ(gm.num_factors(), 7);
    vec_t buffer(2);

    SUBCASE("binary-only"){
        // no unaries
        vec_t fis{4,5,6};

        for(auto fi: fis){
            CHECK_EQ(gm[fi].arity(), 2);
        }

        SUBCASE("all-equal"){
            vec_t labels{0,0,0,0};
            CHECK_EQ(opengm::detail::evaluate_factors(gm,labels, fis.begin(), fis.end(), buffer), doctest::Approx(0.0f));
        }
        SUBCASE("not-all-equal"){
            vec_t labels{0,1,2,3};
            CHECK_GT(opengm::detail::evaluate_factors(gm,labels, fis.begin(), fis.end(), buffer),0);
        }
    }

    SUBCASE("all"){
        vec_t var{0,1,2,3};
        vec_t fis{0,1,2,3,4,5,6};
        vec_t labels(4);
        gm.space().for_each_state(var.begin(), var.end(), labels, [&](auto && labels){

            auto ef = opengm::detail::evaluate_factors(gm,labels, fis.begin(), fis.end(), buffer);
            auto egm = gm.evaluate(labels.begin(), labels.end());
            CHECK_EQ(ef, doctest::Approx(egm));

        });
    }


}



TEST_SUITE_END(); // end of testsuite gm
