#pragma once


template<class gm_type>
auto to_binary(const gm_type & gm)
{
    // how many variables do we need to encode binary
    std::vector<std::size_t> num_var_to_encode_binary(gm.num_variables());
    for(auto vi=0; vi<gm.num_variables(); ++vi)
    {
        num_var_to_encode_binary[vi] = static_cast<std::size_t>(std::ceil(std::log2(gm.num_labels(vi))));
    }
    auto num_binary_var = std::accumulate(num_var_to_encode_binary.begin(), num_var_to_encode_binary.end(), 0);

    // space for the binary graphical model
    UniformSpace<2> binary_space(num_binary_var);


    for(auto && factor : gm)
    {
        auto && variables = factor.variables();
        const auto arity = factor.arity();

        // arity of the new factor
        auto binary_arity = 0;
        for(auto vi : variables)
        {
            binary_arity += num_var_to_encode_binary[vi];
        }

        // all variables of that factor where binary in the first place
        if(binary_arity == arity)
        {

        }
        // at least one variable had more than 2 variables
        else
        {
            // translate the factor to the binary factor with a bigger arity
        }
    }
}