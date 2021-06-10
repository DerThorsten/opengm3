#include <doctest.h>
#include "utils.hpp"

#include "opengm/to_quadratic.hpp"
#include "opengm/toy_models.hpp"



TEST_SUITE_BEGIN("gm");

TEST_CASE("to_quadratic"){

    auto n_variables = 4;
    auto n_factors = 10;
    auto min_num_labels =  2;
    auto max_num_labels = min_num_labels;
    auto seed = 42;
    auto min_arity = 1;
    auto max_arity = 4;
    auto gm = opengm::RandomModel<>(n_variables, n_factors,
        min_num_labels, max_num_labels,
        min_arity, max_arity, seed)();
    using gm_type = std::decay_t<decltype(gm)>;
    using label_type = typename gm_type::label_type;
    using vec_t = std::vector<std::size_t>;
    auto quadratic_gm = to_quadratic(gm);

    CHECK_GE(quadratic_gm.num_variables(), gm.num_variables());
    if(gm.max_arity() >= 2)
    {
        CHECK_EQ(quadratic_gm.max_arity(), 2);
    }

}



TEST_SUITE_END(); // end of testsuite gm
