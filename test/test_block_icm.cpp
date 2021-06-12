#include <doctest.h>

#include "utils.hpp"

#include "opengm/toy_models.hpp"


#include "opengm/minimizer/bp.hpp"
#include "opengm/minimizer/block_icm.hpp"

TEST_SUITE_BEGIN("gm");

namespace cond = opengm::condition;







TEST_CASE("BlockIcm"){

    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::BlockIcm<gm_type>;
        typename minimizer_type::settings_type settings;


        // gm type for the fusion move solver
        using sub_gm_type = typename minimizer_type::sub_gm_type;

        settings.minimizer_factory = opengm::make_shared_factory<opengm::BeliefPropergation<sub_gm_type>>([&](auto & settings){
            settings.num_iterations = 20;
            settings.damping = 0.5;
        });
        return minimizer_type(gm, settings);
    };

    {
        auto nx = 20;
        auto ny = 20;
        auto n_labels =  5;

        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels), factory, cond::NoCondition(), 1);
    }
    {
        auto n_variables = 10;
        auto n_factors = 10;
        auto min_num_labels = 2;
        auto max_num_labels = 3;
        auto min_arity = 1;
        auto max_arity = 4;

        opengm::testing(opengm::RandomModel<>(n_variables, n_factors,min_num_labels, max_num_labels,min_arity, max_arity), factory, cond::NoCondition(), 2);
    }
}





TEST_SUITE_END(); // end of testsuite gm
