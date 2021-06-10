#include <doctest.h>

#include "utils.hpp"

#include "opengm/toy_models.hpp"

#include "opengm/minimizer/icm.hpp"
#include "opengm/minimizer/utils/label_fuser.hpp"

TEST_SUITE_BEGIN("gm");

namespace cond = opengm::condition;




TEST_CASE("LabelFuser"){


    std::random_device rd; 
    std::mt19937 gen(42); // these can be global and/or static, depending on how you use random elsewhere



    auto nx = 5;
    auto ny = 5;
    auto n_labels =  4;




    std::uniform_int_distribution<> label_dist(0, n_labels-1);


    auto gm = opengm::RandomPottsGrid(nx, ny, n_labels)();
    using gm_type = std::decay_t<decltype(gm)>;
    using label_fuser_type = opengm::LabelFuser<gm_type>;
    using label_fuser_settings_type = typename label_fuser_type::settings_type;
    using fuse_gm_type = typename label_fuser_type::fuse_gm_type;
    using labels_vector_type = typename gm_type::labels_vector_type;

    using fuse_minimizer_type = opengm::Icm<fuse_gm_type>;
    using fuse_minimizer_factory_type = opengm::MinimizerFactory<fuse_minimizer_type>;

    // the factory
    auto fuse_minimizer_factory = std::make_shared<fuse_minimizer_factory_type>();

    label_fuser_settings_type label_fuser_settings;
    label_fuser_type label_fuser(gm, label_fuser_settings);


    labels_vector_type labels_a(gm.num_variables());
    labels_vector_type labels_b(gm.num_variables());
    labels_vector_type labels_fused(gm.num_variables());

    // random labels
    std::generate(labels_a.begin(), labels_a.end(), [&](){ return label_dist(gen); });
    std::generate(labels_b.begin(), labels_b.end(), [&](){ return label_dist(gen); });

    // fuse
    label_fuser.fuse(labels_a, labels_b, labels_fused);



}   


TEST_SUITE_END(); // end of testsuite gm
