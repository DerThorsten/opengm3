#pragma once

#include <vector>

#include "opengm/opengm_config.hpp"
#include "opengm/graphical_model.hpp"
#include "opengm/graphical_model.hpp"
#include "opengm/minimizer/utils/fix-fusion/higher-order-energy.hpp"
#include "opengm/space.hpp"

namespace opengm::detail
{
    template <typename REAL>
    struct OpengmQuadraticRep {

    
        using value_type = REAL;
        using space = StaticNumLabelsSpace<label_type,2>;
        using gm_type = GraphicalModel<space, value_type>;
        using unary_tensor_type  = BinaryMultilinearTensor<value_type, 1>;
        using binary_tensor_type = BinaryMultilinearTensor<value_type, 2>;

        typedef int NodeId;

        OpengmQuadraticRep()
        :   m_nodes(0),
            m_gm()
        {
        }


        NodeId AddNode(int n = 1)
        {
            m_nodes += n;
            m_gm.space().resize(m_nodes);
            return m_nodes;
        }
        void SetMaxEdgeNum(int )
        {
            // we might want to call reserve for the vector
            // of factors/tensors of the graphical model
        }

        void AddUnaryTerm(NodeId n, REAL E0, REAL E1)
        {
            // we know that all other coefficients are zero!
            m_gm.add_unary_factor(std::make_unique<unary_tensor_type>(E1), n);
        }
        void AddPairwiseTerm(NodeId n1, NodeId n2, REAL E00, REAL E01, REAL E10, REAL E11)
        {
            // we know that all other coefficients are zero!
            m_gm.add_factor(std::make_unique<binary_tensor_type>(E11), {n1, n2});
        }


        int m_nodes;
        gm_type m_gm;
    };
}

namespace opengm
{

    template<class GM, std::size_t max_allowed_arity = 10>
    class ToQuadratic
    {
    public:
        using gm_type = GM;

        using value_type = typename gm_type::value_type;
        using hoe_type =  HigherOrderEnergy<value_type, max_allowed_arity>;
        using quadratic_gm_type = typename detail::OpengmQuadraticRep<value_type>::gm_type;
        auto static to_quadratic(const gm_type & gm)
        {
            if(gm.space().max_num_labels() != 2)
            {
                throw std::runtime_error("to_quadratic only works for models where the max_num_labels == 2");
            }

            const auto max_arity = gm.max_arity();
            if(gm.max_arity() > max_allowed_arity)
            {
                throw std::runtime_error("arity is bigger than max_allowed_arity");
            }
            hoe_type hoe;
            hoe.AddVars(gm.num_variables());


            const unsigned int max_num_assigments = 1 << max_arity;
            std::vector<value_type> coeffs_array(max_num_assigments);
            std::vector<label_type> clique_labels_array(max_arity);

            for(auto && factor : gm)
            {
                auto && variables = factor.variables();
                const auto arity = factor.arity();
                if(arity == 1)
                {
                    hoe.AddUnaryTerm(variables[0], factor(1) - factor(0));
                }
                else if(arity >= 2)
                {
                    unsigned int num_assigments = 1 << arity;
                    auto coeffs = coeffs_array.data();
                    for (unsigned int subset = 1; subset < num_assigments; ++subset)
                    {
                        coeffs[subset] = 0;
                    }

                    // For each boolean assignment, get the clique energy at the
                    // corresponding labeling
                    auto clique_labels = clique_labels_array.data();
                    for (unsigned int assignment = 0; assignment < num_assigments; ++assignment)
                    {
                        for (unsigned int i = 0; i < arity; ++i)
                        {
                            clique_labels[i] = (assignment & (1 << i) ? 0 : 1);
                        }
                        const auto energy = factor[clique_labels];
                        for (unsigned int subset = 1; subset < num_assigments; ++subset)
                        {
                            if (assignment & ~subset)
                            {
                                continue;
                            }
                            else
                            {
                                int parity = 0;
                                for (unsigned int b = 0; b < arity; ++b)
                                {
                                    parity ^= (((assignment ^ subset) & (1 << b)) != 0);
                                }
                                coeffs[subset] += parity ? -energy : energy;
                            }
                        }
                    }
                    typename hoe_type::VarId vars[max_arity];
                    for (unsigned int subset = 1; subset < num_assigments; ++subset)
                    {
                        int degree = 0;
                        for (unsigned int b = 0; b < arity; ++b)
                        {
                            if (subset & (1 << b)) {
                                vars[degree++] = variables[b];
                            }
                        }
                        std::sort(vars, vars + degree);
                        hoe.AddTerm(coeffs[subset], degree, vars);
                    }
                }
            }

            detail::OpengmQuadraticRep<value_type> opengm_quadratic;
            hoe.ToQuadratic(opengm_quadratic);
            return std::move(opengm_quadratic.m_gm);
        }
    };


    template<class gm_type, std::size_t max_allowed_arity = 10>
    auto to_quadratic(const gm_type & gm)
    {
        return ToQuadratic<gm_type, max_allowed_arity>::to_quadratic(gm);
    }
}