#pragma once

#include <queue>
#include <algorithm>

#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/minimizer/minimizer_base.hpp"
#include "opengm/minimizer/utils/conditioned_submodel.hpp"
#include "opengm/from_gm_factory.hpp"

#include "opengm/minimizer/factor_icm.hpp"


namespace opengm{

template<class GM>
class BlockIcm;

template<class GM>
class BlockGeneratorBase
{
public:
    using gm_type = GM;
    using block_icm_type = BlockIcm<GM>;
    using variables_set_type = std::set<std::size_t>;

    virtual bool generate(block_icm_type * block_icm) = 0;

};

template<class gm_type>
using BlockGeneratorFactoryBase = detail::FromGmFactoryBase<BlockGeneratorBase<gm_type>, gm_type>;


template<class block_generator_type>
using BlockGeneratorFactory = detail::FromGmFactory<
    block_generator_type,
    BlockGeneratorBase<typename block_generator_type::gm_type>,
    typename block_generator_type::gm_type,
    typename block_generator_type::settings_type
>;
// template<class GM>
// class BlockGeneratorFactoryBase{
// public:
//     using gm_type = GM;
//     std::unique_ptr<GM>
// };


template<class GM>
class BlockIcm : public MinimizerCrtpBase<GM, BlockIcm<GM> >{
public:
    using gm_type = GM;
    using base_type = MinimizerBase<GM>;
    using value_type = typename GM::value_type;
    using label_type = typename GM::label_type;
    using labels_vector_type = typename base_type::labels_vector_type;
    using minimizer_callback_base_ptr_type = typename base_type::minimizer_callback_base_ptr_type;

    using conditioned_submodel_builder_type = detail::ConditionedSubmodelBuilder<gm_type>;
    using sub_gm_type = typename conditioned_submodel_builder_type::sub_gm_type; 



    using default_sub_minimizer_type = FactorIcm<sub_gm_type>;
    using default_sub_minimier_factory_type = MinimizerFactory<default_sub_minimizer_type>;

    using sub_minimizer_factory_base_type = MinimizerFactoryBase<sub_gm_type>;
    using sub_minimizer_factory_ptr_type = std::shared_ptr<sub_minimizer_factory_base_type>;


    using base_type::minimize;


    struct settings_type : public SolverSettingsBase{
        public:
            sub_minimizer_factory_ptr_type minimizer_factory;
            //value_type eps{1e-7};
    };


    BlockIcm(const GM & gm, const settings_type & settings = settings_type())
    :   m_gm(gm),
        m_labels(gm.num_variables(),0),
        m_energy(m_gm.evaluate(m_labels))
    {
    }

    std::string name() const override{
        return "BlockIcm";
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
        return true;
    }
    void set_starting_point(const label_type  * labels)override{
        m_labels.assign(labels, labels + m_gm.num_variables());
    }

    void minimize(minimizer_callback_base_ptr_type minimizer_callback_base_ptr)override{
        auto callback = callback_wrapper(this, minimizer_callback_base_ptr);

        m_energy = m_gm.evaluate(m_labels);

    }
private:


    const GM & m_gm;
    std::vector<label_type> m_labels;
    value_type m_energy;

};

template<class GM>
using BlockIcmFactory = MinimizerFactory<BlockIcm<GM>>;
}