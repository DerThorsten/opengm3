#include <doctest.h>

#include "utils.hpp"
#include "opengm/tensors.hpp"



TEST_SUITE_BEGIN("gm");

TEST_CASE("XArrayTensor"){

    using tensor_type = opengm::XArrayTensor<int>;
    using xarray_shape = typename tensor_type::xshape_type;
    tensor_type tensor(xarray_shape({2,2}));
}



TEST_CASE("TestBind"){

    std::size_t fixed_pos[1];
    std::size_t fixed_labels[1];



    opengm::Potts2Tensor<int> tensor(3, 1);
    fixed_pos[0] = 0;
    fixed_labels[0] = 2;
    auto binded_tensor = tensor.bind(
        gsl::span<const std::size_t>(fixed_pos, 1),
        gsl::span<const std::size_t>(fixed_labels, 1)
    );
    CHECK_EQ(tensor(1,2), 1);

    CHECK_EQ(binded_tensor->arity() , 1);
    CHECK_EQ(binded_tensor->shape(0) , tensor.shape(1));

    std::size_t l=0;
    CHECK_EQ(binded_tensor->operator[](&l), tensor(2,0));
    l=1;
    CHECK_EQ(binded_tensor->operator[](&l), tensor(2,1));
    l=2;
    CHECK_EQ(binded_tensor->operator[](&l), tensor(2,2));


}







TEST_SUITE_END(); // end of testsuite gm
