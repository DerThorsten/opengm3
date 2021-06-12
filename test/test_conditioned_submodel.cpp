#include <doctest.h>

#include "utils.hpp"
#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/toy_models.hpp"
#include "opengm/minimizer/brute_force_naive.hpp"



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

    conditioned_submodel_builder.condition({0}, labels, [&](auto && sub_gm){
        CHECK_EQ(sub_gm.num_variables(), 1);
        for(auto && factor : sub_gm)
        {
            if(factor.arity() == 1)
            {
                auto a0 = factor(0);
                auto a1 = factor(1);
                auto a2 = factor(2);
            }
            if(factor.arity() == 2)
            {
                auto a00 = factor(0,0);
                auto a11 = factor(1,2);
                auto a22 = factor(2,2);
            }
        }
    });
}


TEST_CASE("ConditionedSubmodel2"){



    
    auto n_variables = 10;
    auto n_factors = 10;
    auto min_num_labels = 2;
    auto max_num_labels = 3;
    auto min_arity = 1;
    auto max_arity = 4;

    auto gen = opengm::RandomModel<>(n_variables, n_factors,min_num_labels, max_num_labels,min_arity, max_arity);
    

    for(auto i=0;i<10;++i)
    {

        auto gm = gen();

        using gm_type = std::decay_t<decltype(gm)>;
        using labels_vector_type = typename gm_type::labels_vector_type;
        auto conditioned_submodel_builder = opengm::detail::conditioned_submodel_builder(gm);

        labels_vector_type labels(gm.num_variables(), 1);

        conditioned_submodel_builder.condition({0}, labels, [&](auto && sub_gm){
            using sub_gm_type = std::decay_t<decltype(sub_gm)>;
            CHECK_EQ(sub_gm.num_variables(), 1);
            for(auto && factor : sub_gm)
            {
                if(factor.arity() == 1)
                {
                    auto a0 = factor(0);
                    auto a1 = factor(1);
                    auto a2 = factor(2);
                }
                if(factor.arity() == 2)
                {
                    auto a00 = factor(0,0);
                    auto a11 = factor(1,2);
                    auto a22 = factor(2,2);
                }
                if(factor.arity() == 3)
                {
                    auto a00 = factor(0,0,0);
                    auto a11 = factor(1,2,0);
                    auto a22 = factor(2,2,0);
                }
            }


            auto minimizer_factory = opengm::make_shared_factory<opengm::BruteForceNaive<sub_gm_type>>([&](auto & settings){
            });
            auto minimizer = minimizer_factory->create(gm);
            minimizer->minimize();
        });
    }
}



TEST_SUITE_END(); // end of testsuite gm
