#include <doctest.h>

#include "utils.hpp"

#include "opengm/toy_models.hpp"


#include "opengm/minimizer/icm.hpp"
#include "opengm/minimizer/factor_icm.hpp"
#include "opengm/minimizer/brute_force_naive.hpp"
#include "opengm/minimizer/bp.hpp"
#include "opengm/minimizer/dynamic_programming.hpp"
#include "opengm/minimizer/self_fusion.hpp"
#include "opengm/minimizer/block_icm.hpp"

TEST_SUITE_BEGIN("gm");

namespace cond = opengm::condition;



TEST_CASE("SelfFusion"){

    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::SelfFusion<gm_type>;
        typename minimizer_type::settings_type settings;


        // gm type for the fusion move solver
        using sub_gm_type = typename minimizer_type::sub_gm_type;

        settings.minimizer_factory = opengm::make_shared_factory<opengm::BeliefPropergation<gm_type>>([&](auto & settings){
            settings.num_iterations = 20;
            settings.damping = 0.5;
        });
        settings.fuse_minimizer_factory = opengm::make_shared_factory<opengm::FactorIcm<sub_gm_type>>();
        return minimizer_type(gm, settings);
    };

    {
        auto nx = 20;
        auto ny = 20;
        auto n_labels =  5;

        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels), factory, cond::NoCondition(), 1);
    }
}

TEST_CASE("DynamicProgramming"){

    // lambda as generic factory
    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::DynamicProgramming<gm_type>;
        typename minimizer_type::settings_type settings;
        return minimizer_type(gm, settings);
    };
    {
        auto n_variables = 9;
        auto n_labels =  2;
        opengm::testing(opengm::RandomPottsChain(n_variables, n_labels), factory, cond::Optimal(), 1000);
    }
    {
        auto n_variables = 6;
        auto n_labels =  3;
        opengm::testing(opengm::RandomPottsChain(n_variables, n_labels), factory, cond::Optimal(), 1000);
    }
    {
        auto nx = 3;
        auto ny = 2;
        auto n_labels =  3;
        opengm::testing(opengm::RandomPottsGridFan(nx, ny, n_labels), factory, cond::Optimal(), 1000);
    }
    {
        auto nx = 3;
        auto ny = 3;
        auto n_labels =  2;
        opengm::testing(opengm::RandomPottsGridFan(nx, ny, n_labels), factory, cond::Optimal(), 1000);
    }
}



TEST_CASE("BeliefPropergation"){

    // lambda as generic factory
    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::BeliefPropergation<gm_type>;
        typename minimizer_type::settings_type settings;
        settings.damping = 0.5;
        settings.num_iterations = 100;
        return minimizer_type(gm, settings);
    };
    {
        auto nx = 10;
        auto ny = 10;
        auto n_labels =  4;

        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels), factory, cond::NoCondition(), 100);
    }
    {
        auto n_variables = 6;
        auto n_labels =  4;
        opengm::testing(opengm::RandomPottsChain(n_variables, n_labels), factory, cond::NoCondition(), 10);
    }
}

TEST_CASE("Icm"){

    // lambda as generic factory
    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::Icm<gm_type>;
        typename minimizer_type::settings_type settings;
        return minimizer_type(gm, settings);
    };
    {
        auto nx = 10;
        auto ny = 10;
        auto n_labels =  4;

        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels), factory, cond::KOptimal<1>(), 100);
    }
    {
        auto n_variables = 6;
        auto n_labels =  4;
        opengm::testing(opengm::RandomPottsChain(n_variables, n_labels), factory, cond::KOptimal<1>(), 10);
    }
}


TEST_CASE("FactorIcm"){

    // lambda as generic factory
    auto factory =  [](auto && gm){
        using gm_type = std::decay_t<decltype(gm)>;
        using minimizer_type = opengm::FactorIcm<gm_type>;
        typename minimizer_type::settings_type settings;
        return minimizer_type(gm, settings);
    };
    {
        auto nx = 5;
        auto ny = 5;
        auto n_labels =  4;

        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels), factory, cond::KOptimal<2>(), 100);
    }
    {
        auto nx = 2;
        auto ny = 1;
        auto n_labels =  50;

        opengm::testing(opengm::RandomPottsChain(nx, ny, n_labels), factory, cond::Optimal(), 100);
    }
    {
        auto nx = 2;
        auto ny = 2;
        auto n_labels =  2;
        auto seed = 46347;
        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels, seed), factory, cond::FactorLocalOptimal(), 1);
    }
    {
        auto nx = 5;
        auto ny = 5;
        auto n_labels =  3;
        auto seed = 0;
        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels, seed), factory,  cond::KOptimal<2>(), 100);
    }
    {
        auto nx = 10;
        auto ny = 10;
        auto n_labels =  2;
        auto seed = 0;
        opengm::testing(opengm::RandomPottsGrid(nx, ny, n_labels, seed), factory, cond::FactorLocalOptimal(), 100);
    }
    {
        auto n_variables = 7;
        auto n_labels =  4;
        opengm::testing(opengm::RandomPottsChain(n_variables, n_labels), factory,cond::FactorLocalOptimal(), 100);
    }
}

TEST_SUITE_END(); // end of testsuite gm
