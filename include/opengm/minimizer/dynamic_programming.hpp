#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/factors_of_variables.hpp"
#include "opengm/factors_of_variables.hpp"

namespace opengm{


template<class GM>
class DynamicProgramming : public MinimizerCrtpBase<GM, DynamicProgramming<GM> >{
public:
    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;

    using base_type::minimize;

    struct settings_type : public SolverSettingsBase{
        std::vector<std::size_t> roots;
    };


    DynamicProgramming(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_settings(settings),
        m_factors_of_variables(gm),
        m_energy(),
        m_labels(gm.num_variables(),0),
        m_value_buffer(),
        m_state_buffer(),
        m_value_buffers(gm.num_variables()),
        m_state_buffers(gm.num_variables()),
        m_node_order(gm.num_variables(), std::numeric_limits<std::size_t>::max() ),
        m_ordered_nodes(gm.num_variables(), std::numeric_limits<std::size_t>::max() )
    {

        if(m_gm.max_arity() > 2)
        {
            throw std::runtime_error("This implementation of DynamicProgramming does only support second order models");
        }

        // node order
        std::vector<std::size_t> num_children(m_gm.num_variables(), 0);
        std::vector<std::size_t> node_list;

        std::size_t order_count = 0;
        std::size_t var_count = 0;
        std::size_t root_count = 0;


        constexpr auto mxval = std::numeric_limits<std::size_t>::max() ;
        while(var_count < m_gm.num_variables() && order_count < m_gm.num_variables())
        {
            if(root_count<m_settings.roots.size())
            {
                m_node_order[m_settings.roots[root_count]] = order_count++;
                node_list.push_back(m_settings.roots[root_count]);
                ++root_count;
            }
            else if(m_node_order[var_count]==std::numeric_limits<std::size_t>::max())
            {
                m_node_order[var_count] = order_count++;
                node_list.push_back(var_count);
            }
            ++var_count;
            while(node_list.size()>0){
                size_t node = node_list.back();
                node_list.pop_back();
                for(auto && fid: m_factors_of_variables[node])
                {
                    auto && factor = m_gm[fid];
                    auto && variables = factor.variables();
                    if( factor.arity() == 2 ){
                        if(variables[1] == node && m_node_order[variables[0]]==mxval ){
                            m_node_order[variables[0]] = order_count++;
                            node_list.push_back(variables[0]);
                            ++num_children[node];
                        }
                        if( variables[0] == node && m_node_order[variables[1]]==mxval ){
                            m_node_order[variables[1]] = order_count++;
                            node_list.push_back(variables[1]);
                            ++num_children[node];
                        }
                    }
                }
            }
        }

        auto buffer_size_values = 0;
        auto buffer_size_states = 0;
        for(std::size_t vi=0; vi<m_gm.num_variables();++vi)
        {
            buffer_size_values += m_gm.num_labels(vi);
            buffer_size_states += gm.num_labels(vi) * num_children[vi];
        }
        m_value_buffer.resize(buffer_size_values);
        m_state_buffer.resize(buffer_size_states);

        auto value_ptr =  m_value_buffer.data();
        auto state_ptr =  m_state_buffer.data();

        for(std::size_t vi=0; vi<m_gm.num_variables();++vi)
        {
            m_value_buffers[vi] = value_ptr;
            value_ptr += gm.num_labels(vi);
            m_state_buffers[vi] = state_ptr;
            state_ptr +=  gm.num_labels(vi) * num_children[vi];
            m_ordered_nodes[m_node_order[vi]] = vi;
        }
    }



    std::string name() const override{
        return "DynamicProgramming";
    }

    const gm_type & gm() const override{
        return m_gm;
    }

    const labels_vector_type & best_labels()override{
        return m_labels;
    }
    const labels_vector_type & current_labels() override{
        return m_labels;
    }

    value_type best_energy()  override{
        return m_energy;
    }
    value_type current_energy() override {
        return m_energy;
    }


    bool can_start_from_starting_point() override{
        return false;
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override {
        m_energy = m_gm.evaluate(m_labels);
        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);


        for (std::size_t i = 1; i <= m_gm.num_variables(); ++i)
        {
            // std::cout<<" ii "<< i<<"\n";
            const auto node = m_ordered_nodes[m_gm.num_variables() - i];

            std::fill(m_value_buffers[node], m_value_buffers[node] + m_gm.num_labels(node), value_type(0));

            // accumulate messages
            std::size_t children_counter = 0;
            for(auto fid: m_factors_of_variables[node])
            {
                auto && factor = m_gm[fid];

                // unary
                if (factor.arity() == 1)
                {
                    factor.add_values(m_value_buffers[node]);
                }

                //pairwise
                if (factor.arity() == 2)
                {

                    auto && vars =  factor.variables();
                    if (vars[0] == node && m_node_order[vars[1]] > m_node_order[node])
                    {
                        const auto node2 = vars[1];
                        label_type s;
                        value_type v;
                        for(auto l0=0; l0<m_gm.num_labels(node); ++l0)
                        {
                            v=std::numeric_limits<value_type>::infinity();
                            for(auto l1=0; l1<m_gm.num_labels(node2); ++l1)
                            {
                                const auto factor_value = factor(l0, l1);
                                const auto v2 = factor_value + m_value_buffers[node2][l1];
                                if(v2 < v)
                                {
                                    v = v2;
                                    s = l1;
                                }
                            }
                            m_state_buffers[node][children_counter * m_gm.num_labels(node) + l0] = s;
                            m_value_buffers[node][l0] += v;
                        }
                        ++children_counter;

                    }
                    if (vars[1] == node && m_node_order[vars[0]] > m_node_order[node])
                    {
                        const auto node2 = vars[0];
                        label_type s;
                        value_type v;
                        for (auto l1 = 0; l1 < m_gm.num_labels(node); ++l1) {
                            v=std::numeric_limits<value_type>::infinity();
                            for (auto l0 = 0; l0 < m_gm.num_labels(node2); ++l0) {
                                const auto factor_value = factor(l0, l1);
                                const auto v2 = factor_value + m_value_buffers[node2][l0];
                                if (v2 < v) {
                                    v = v2;
                                    s = l0;
                                }
                            }
                            m_state_buffers[node][children_counter * m_gm.num_labels(node) + l1] = s;
                            m_value_buffers[node][l1] += v;
                        }
                        ++children_counter;
                    }
                }
            }
        }
        this->compute_labels();
        m_energy = m_gm.evaluate(m_labels);
        callback();
    }
private:

    void compute_labels(){

        std::vector<std::size_t> node_list;
        std::fill(m_labels.begin(), m_labels.end(), std::numeric_limits<label_type>::max());
        std::size_t var = 0;
        while (var < m_gm.num_variables())
        {
            if (m_labels[var] == std::numeric_limits<label_type>::max())
            {
                value_type v = std::numeric_limits<value_type>::infinity();
                for (std::size_t i = 0; i < m_gm.num_labels(var); ++i) {
                    if(m_value_buffers[var][i] < v){
                        v = m_value_buffers[var][i];
                        m_labels[var] = i;
                    }
                }
                node_list.push_back(var);
            }
            ++var;
            while (node_list.size() > 0)
            {
                std::size_t node = node_list.back();
                std::size_t children_counter = 0;
                node_list.pop_back();


                for(auto && fid: m_factors_of_variables[node])
                {
                    auto && factor = m_gm[fid];
                    auto && vars = factor.variables();

                    if (factor.arity() == 2)
                    {
                        if (vars[1] == node && m_node_order[vars[0]] > m_node_order[node]) {
                            m_labels[vars[0]] = m_state_buffers[node][children_counter * m_gm.num_labels(node) + m_labels[node]];
                            node_list.push_back(vars[0]);
                            ++children_counter;
                        }
                        if (vars[0] == node && m_node_order[vars[1]] > m_node_order[node]) {
                            m_labels[vars[1]] = m_state_buffers[node][children_counter * m_gm.num_labels(node) + m_labels[node]];
                            node_list.push_back(vars[1]);
                            ++children_counter;
                        }
                    }
                }
            }
        }
    }

    const GM & m_gm;
    settings_type m_settings;
    FactorsOfVariables<GM> m_factors_of_variables;
    value_type m_energy;
    std::vector<label_type> m_labels;


    std::vector<value_type> m_value_buffer;
    std::vector<label_type> m_state_buffer;
    std::vector<value_type*> m_value_buffers;
    std::vector<label_type*> m_state_buffers;
    std::vector<std::size_t> m_node_order;
    std::vector<std::size_t> m_ordered_nodes;
};

template<class GM>
using DynamicProgrammingFactory = MinimizerFactory<DynamicProgramming<GM>>;

}