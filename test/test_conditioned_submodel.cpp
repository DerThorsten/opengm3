#include <doctest.h>

#include "utils.hpp"
#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/toy_models.hpp"



TEST_SUITE_BEGIN("gm");




TEST_CASE("TreeFinder"){


    auto gm = opengm::RandomPottsGrid(10/*nx*/,10/*ny*/,3/*n_labels*/)();
    using gm_type = std::decay_t<decltype(gm)>;
    using labels_vector_type = typename gm_type::labels_vector_type;
    auto conditioned_submodel_builder = opengm::detail::conditioned_submodel_builder(gm);

    using tree_finder_type = opengm::detail::TreeFinder<gm_type>;
    using factors_of_variables_type = typename tree_finder_type::factors_of_variables_type;

    factors_of_variables_type factors_of_variables(gm);
    tree_finder_type tree_finder(gm, factors_of_variables);

    tree_finder.find(0);
}



TEST_CASE("ConditionedSubmodel"){


    auto gm = opengm::RandomPottsGrid(10/*nx*/,10/*ny*/,3/*n_labels*/)();
    using gm_type = std::decay_t<decltype(gm)>;
    using labels_vector_type = typename gm_type::labels_vector_type;
    auto conditioned_submodel_builder = opengm::detail::conditioned_submodel_builder(gm);

    labels_vector_type labels(gm.num_variables(), 1);

    conditioned_submodel_builder.condition({0,1,2}, labels, [&](auto && sub_gm){
        CHECK_EQ(sub_gm.num_variables(), 3);
    });
}



TEST_SUITE_END(); // end of testsuite gm
