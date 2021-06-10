#pragma once

#include "opengm/opengm_config.hpp"
#include "opengm/graphical_model.hpp"
#include "opengm/graphical_model.hpp"
#include "opengm/space.hpp"

namespace opengm::detail
{
    template <typename REAL>
    class DummyQuadraticRep {

        private:
            using value_type = REAL;
            using space = StaticNumLabelsSpace<label_type,2>;
            using gm_type = GraphicalModel<space, value_type>;
            using unary_tensor_type  = BinaryMultilinearTensor<value_type, 1>;
            using binary_tensor_type = BinaryMultilinearTensor<value_type, 2>;

        public:
            typedef int NodeId;

            DummyQuadraticRep()
            :   m_nodes(0),
                m_gm()
            {
            }


            NodeId AddNode(int n)
            {
                m_nodes += n;
                m_gm.space().resize(m_nodes);
                return m_nodes;
            }
            void SetMaxEdgeNum(int n)
            {
            }

            // void AddConstantTerm(REAL c)
            // {
            //     throw std::runtime_error("dead code reached");
            // }
            // void AddUnaryTerm(NodeId n, REAL coeff)
            // {
            //     throw std::runtime_error("dead code reached");
            // }
            void AddUnaryTerm(NodeId n, REAL E0, REAL E1)
            {
                m_gm.add_unary_factor(std::make_unique<unary_tensor_type>(E1), n);
            }
            // void AddPairwiseTerm(NodeId n1, NodeId n2, REAL coeff)
            // {
            //     throw std::runtime_error("dead code reached");
            // }
            void AddPairwiseTerm(NodeId n1, NodeId n2, REAL E00, REAL E01, REAL E10, REAL E11)
            {
                // we know that all other coefficients are zero!
                m_gm.add_factor(std::make_unique<binary_tensor_type>(E11), {n1, n2});
            }

            // void MergeParallelEdges() {
            //     throw std::runtime_error("dead code reached");
            // }
            // void Solve() {
            //     auto minimizer =  m_minimizer_factory.create(m_gm);
            //     minimizer->optimize();
            // }

            int GetNodeNum() { return m_nodes; }
            int GetLabel(NodeId n) { return -1; }

        protected:
            int m_nodes;
            gm_type m_gm;
            //minimizer_factory_ptr_type m_minimizer_factory;

    };
}