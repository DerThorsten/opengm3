#pragma once



#include <array>
#include <vector>
#include <algorithm>

#include "opengm/opengm_config.hpp"
#include "opengm/utils.hpp"
#include "opengm/crtp_base.hpp"

namespace opengm {


    template<class derived>
    class SpaceTraits{

    };

    template<class derived>
    class SpaceBase : public CrtpBase<derived>{
    public:
        using derived_t = derived;
        using space_traits_t = SpaceTraits<derived_t>;
        using label_type = typename space_traits_t::label_type;
        using subspace_type = typename space_traits_t::subspace_type;

        auto max_num_labels()const{
            label_type mx = 0;
            for(auto i=0; i<this->derived_cast().size(); ++i){
                mx = std::max(mx, this->derived_cast()[i]);
            }
            return mx;
        }
        auto is_uniform_space()const{
            label_type l =  this->derived_cast()[0];
            for(auto i=1; i<this->derived_cast().size(); ++i){
                if(this->derived_cast()[i] != l){
                    return false;
                }
            }
            return true;
        }


        template<class LABELS, class F>
        void for_each_state(
            LABELS & labels,
            F && f
        ) const
        {
            detail::for_each_state(this->derived_cast().size(),
                this->derived_cast(), labels, std::forward<F>(f));
        }


        template<class ITER, class LABELS, class F>
        void for_each_state(
            ITER vars_begin,
            ITER vars_end,
            LABELS & labels,
            F && f
        ) const
        {
            // initialize labels with zero
            std::for_each(vars_begin, vars_end, [&](auto var){labels[var] = 0;});

            const auto n_var = std::distance(vars_begin, vars_end);

            for (;;)
            {

                // do smth with state
                f(labels);

                // increment buffered state
                for (size_t j = n_var-1; j >= 0; --j)
                {
                    const auto vi = vars_begin[j];
                    if (labels[vi] < this->derived_cast()[vi] - 1)
                    {
                        ++labels[vi];
                        break;
                    }
                    else
                    {
                        if(j != 0)
                        {
                            labels[vi] = 0;
                        }
                        else
                        {
                            goto overflow;
                        }
                    }
                }

            }
            overflow:
            ;
        }
    };



    template<class label_type=label_type>
    class UniformSpace : public SpaceBase<UniformSpace<label_type>>{
    public:
        using self_type =  UniformSpace<label_type>;
        using subspace_type = self_type;

        UniformSpace(const std::size_t n_var = 0, const label_type n_labels = label_type(0))
        :   m_n_var(n_var),
            m_n_labels(n_labels){
        }
        auto operator[](const std::size_t)const{
            return m_n_labels;
        }
        auto size()const{
            return m_n_var;
        }
        auto max_num_labels()const{
            return m_n_labels;
        }
        constexpr auto is_uniform_space()const{
            return true;
        }

        template<class VI_ITER>
        self_type subspace(VI_ITER begin, VI_ITER end)const{
            return self_type(std::distance(begin, end), m_n_labels);
        }

    private:
        std::size_t m_n_var;
        label_type m_n_labels;
    };

    template<class label_type_t>
    class SpaceTraits< UniformSpace<label_type_t>>{
    public:
        using label_type = label_type_t;
        using subspace_type = UniformSpace<label_type_t>;
    };





    template<class label_type>
    class ExplicitSpace : public SpaceBase<ExplicitSpace<label_type>>{
    public:
        using self_type =  ExplicitSpace<label_type>;
        using subspace_type = self_type;

        ExplicitSpace(std::size_t num_var = 0, label_type num_labels = 0)
        :   m_space(num_var, num_labels)
        {
        }

        template<class T>
        ExplicitSpace(std::initializer_list<T> values)
        :   ExplicitSpace(values.begin(), values.end())
        {
        }

        template<class ITER>
        ExplicitSpace(ITER val_begin, ITER val_end)
        :   m_space(val_begin, val_end)
        {
        }

        auto operator[](const std::size_t var)const{
            return m_space[var];
        }
        auto size()const{
            return m_space.size();
        }

        template<class VI_ITER>
        self_type subspace(VI_ITER begin, VI_ITER end)const{
            self_type ret(std::distance(begin, end));
            auto i=0;
            while(begin != end){
                ret.m_space[i] = m_space[*begin];
                ++begin;
            }
            return ret;
        }

    private:
        std::vector<std::size_t> m_space;
    };

    template<class label_type_t>
    class SpaceTraits< ExplicitSpace<label_type_t>>{
    public:
        using label_type = label_type_t;
        using subspace_type = ExplicitSpace<label_type_t>;
    };



    template<class label_type,label_type num_labels>
    class StaticNumLabelsSpace : public SpaceBase<StaticNumLabelsSpace<label_type, num_labels>>{
    public:
        using self_type = StaticNumLabelsSpace<label_type, num_labels>;
        using subspace_type = self_type;
        StaticNumLabelsSpace(const std::size_t n_var = 0)
        :   m_n_var(n_var){
        }

        auto operator[](const std::size_t var)const{
            return num_labels;
        }
        auto size()const{
            return m_n_var;
        }
        constexpr auto max_num_labels()const{
            return num_labels;
        }
        constexpr auto is_uniform_space()const{
            return true;
        }
        void resize(const std::size_t n_var){
            m_n_var = n_var;
        }

        template<class VI_ITER>
        self_type subspace(VI_ITER begin, VI_ITER end)const{
            return self_type(std::distance(begin, end));
        }

    private:
        std::size_t m_n_var;
    };


    template<class label_type_t, label_type_t num_labels>
    class SpaceTraits< StaticNumLabelsSpace<label_type_t, num_labels>>{
    public:
        using label_type = label_type_t;
        using subspace_type = StaticNumLabelsSpace<label_type_t, num_labels>;
    };

    template<class label_type, std::size_t num_var, label_type num_labels>
    class StaticSpace : public SpaceBase<StaticSpace<label_type, num_var, num_labels>>{
    public:
        using subspace_type = StaticNumLabelsSpace<label_type, num_labels>;
        auto operator[](const std::size_t var)const{
            return num_labels;
        }
        auto size()const{
            return num_var;
        }
        constexpr auto max_num_labels()const{
            return num_labels;
        }
        constexpr auto is_uniform_space()const{
            return true;
        }
        template<class VI_ITER>
        auto subspace(VI_ITER begin, VI_ITER end)const{
            return StaticNumLabelsSpace<label_type, num_labels>(std::distance(begin, end));
        }
    };

    template<class label_type_t, std::size_t num_var, label_type_t num_labels>
    class SpaceTraits< StaticSpace<label_type_t, num_var, num_labels>>{
    public:
        using label_type = label_type_t;
        using subspace_type = StaticNumLabelsSpace<label_type, num_labels>;
    };
}