#pragma once

namespace opengm::detail{

    template<class BASE_TYPE, class GM>
    class FromGmFactoryBase{
    public:
        using gm_type = GM;
        std::unique_ptr<BASE_TYPE> create(const GM & gm) const = 0;
    };

    template<class DERIVED, class BASE_TYPE, class GM, class SETTINGS>
    class FromGmFactory : public FromGmFactoryBase<BASE_TYPE,GM>{
    public:
        using settings_type = SETTINGS;
        FromGmFactory(const settings_type & settings)
        :   m_settings(settings)
        {

        };

        std::unique_ptr<BASE_TYPE> create(const GM & gm) const override{
            return std::make_unique<DERIVED>(gm, m_settings);
        }
    private:
        settings_type m_settings;
    };

}