#pragma once

#include <queue>
#include <algorithm>
#include <map>

#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/minimizer/utils/movemaker.hpp"

namespace opengm{


template<class GM>
class FactorIcm :public MinimizerCrtpBase<GM, FactorIcm<GM> >{
public:
    using gm_type = GM;

    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;
    using movemaker_type = Movemaker<GM>;

    using base_type::minimize;

    struct settings_type : public SolverSettingsBase{
    };

    FactorIcm(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_movemaker(gm),
        m_dirty(),
        m_in_queue(gm.num_factors())
    {

    }

    const gm_type & gm() const override{
        return m_gm;
    }

    const labels_vector_type & best_labels()override{
        return m_movemaker.labels();
    }
    const labels_vector_type & current_labels() override{
        return m_movemaker.labels();
    }
    value_type best_energy()  override{
        return m_movemaker.energy();
    }
    value_type current_energy() override {
        return m_movemaker.energy();
    }
    std::string name() const override{
        return "FactorIcm";
    }

    bool can_start_from_starting_point() override{
        return true;
    }

    void set_starting_point(const label_type  * labels)override{
        m_movemaker.set_labels(labels, labels + m_gm.num_variables());
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override{

        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);


        // start with all factors on the queue
        this->add_all_to_queue();

        for(;;)
        {
            // get next factor
            auto [queue_empty, fi] = this->try_pop_from_queue();
            if(queue_empty)
            {
                break;
            }
            auto && factor = m_gm[fi];



            std::stringstream conf_pre;
            for(auto i=0; i<m_gm.num_variables();++i){
                conf_pre<<this->best_labels()[i]<<" ";
            }


            // move the variables of that factor optimal
            auto && variables = factor.variables();
            this->m_movemaker.move_brute_force_optimal(variables.begin(), variables.end(),[&](auto & involved_factors){

                // this is called *iff* there are changes
                std::set<std::size_t> vars(variables.begin(), variables.end());

                // inform callback since labels changed
                callback();

                for(auto involved_fi : involved_factors)
                {
                    for(auto vi : m_gm[involved_fi].variables())
                    {
                        for(auto other_fi : m_movemaker.factors_of_variables()[vi])
                        {
                            // do *not* add the factor we just optimized
                            // and *not* add factors already in the queue
                            if(other_fi != fi && !m_in_queue[fi])
                            {
                                // if all variables of "other_fi"
                                // are in *variables*, we do not need to
                                // add this factor
                                bool found_all = true;
                                for(auto var : m_gm[other_fi].variables())
                                {
                                    if(vars.find(var) == vars.end())
                                    {
                                        found_all = false;
                                        break;
                                    }
                                }

                                if(!found_all)
                                {
                                    this->try_add_to_queue(other_fi);
                                }
                            }
                        }
                    }
                }
            });
        }
    }
private:

    void add_all_to_queue(){
        for(std::size_t fi=0; fi<m_gm.num_factors(); ++fi){
            m_dirty[m_gm[fi].arity()].push(fi);
            m_in_queue[fi] = true;
        }
    }


    auto try_pop_from_queue(){
        for(auto & kv : m_dirty){
            if(!kv.second.empty()){
                auto fi = kv.second.front();
                m_in_queue[fi] = false;
                kv.second.pop();
                return std::make_pair(false, fi);
            }
        }
        return std::make_pair(true, std::size_t());
    }


    auto try_add_to_queue(const std::size_t fi){
        if(!m_in_queue[fi]){
            m_in_queue[fi] = true;
            m_dirty[m_gm[fi].arity()].push(fi);
        }
    }

    const GM & m_gm;
    movemaker_type m_movemaker;
    std::vector<bool> m_in_queue;
    // arity is the key,    queue the value
    std::map<std::size_t,  std::queue<std::size_t> > m_dirty;
};

template<class GM>
using FactorIcmFactory = MinimizerFactory<FactorIcm<GM>>;

}