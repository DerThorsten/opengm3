#pragma once

#include <vector>
#include <numeric>
#include <cmath>
#include <mutex>

#include "opengm/meta.hpp"
#include "opengm/crtp_base.hpp"
#include "opengm/arity_vector.hpp"
#include "opengm/utils.hpp"
#include "opengm/opengm_config.hpp"

#include <xtensor/xarray.hpp>
#include <gsl-lite/gsl-lite.hpp>

namespace opengm{


namespace detail{

    template<class T>
    inline std::pair<label_type, label_type> arg_2_min(const T * begin, const T * end){
        std::pair<label_type, label_type> res;
        T min0=std::numeric_limits<T>::infinity();
        T min1=min0;
        const auto dist = std::distance(begin, end);
        for(auto i=0; i<dist; ++i){
            const T val = begin[i];
            if(val<min1){
                if(val<min0){

                    min1 = min0;
                    res.second = res.first;
                    min0 = val;
                    res.first = i;
                    continue;
                }
                else{
                    min1 = val;
                    res.second = i;
                }
            }
        }
        return res;
    }

    template<class value_type>
    void potts2_factor_to_variable_messages(
        const label_type num_labels,
        const value_type beta,
        const value_type ** in_messages,
        value_type ** out_messages
    ){
        if(beta>=0){
            const auto min_in_0_beta = *std::min_element(in_messages[0], in_messages[0] + num_labels) + beta;
            const auto min_in_1_beta = *std::min_element(in_messages[1], in_messages[1] + num_labels) + beta;
            for(auto l=0; l<num_labels; ++l){
                out_messages[0][l] = std::min(in_messages[1][l], min_in_1_beta);
                out_messages[1][l] = std::min(in_messages[0][l], min_in_0_beta);
            }
        }
        else{
            auto [a_min_0, a_s_min_0] = detail::arg_2_min(in_messages[0], in_messages[0] + num_labels);
            auto [a_min_1, a_s_min_1] = detail::arg_2_min(in_messages[1], in_messages[1] + num_labels);
            const auto min0 = in_messages[0][a_min_0];
            const auto min1 = in_messages[1][a_min_1];
            const auto min_0_beta = min0 + beta;
            const auto min_1_beta = min1 + beta;
            const auto smin_0_beta = std::min(in_messages[0][a_s_min_0] + beta, min0);
            const auto smin_1_beta = std::min(in_messages[1][a_s_min_1] + beta, min1);
            for(auto l=0; l<num_labels; ++l){
                out_messages[0][l] = a_min_1 != l ?  min_1_beta  : smin_1_beta;
                out_messages[1][l] = a_min_0 != l ?  min_0_beta  : smin_0_beta;
            }
        }
    }

    template<class VT>
    inline void generic_second_order_factor_to_variable_messages(
        const VT * vt,
        const label_type nl_0,
        const label_type nl_1,
        const typename VT::value_type ** in_messages,
        typename VT::value_type ** out_messages
    ){
        using value_type = typename VT::value_type;
        // initialize
        for(auto a=0; a<2; ++a)
        {
            std::fill(out_messages[a], out_messages[a] + (a==0? nl_0 : nl_1), std::numeric_limits<value_type>::infinity());
        }
        // minimize
        for(auto l0=0; l0 < nl_0; ++l0)
        {
            for(auto l1=0; l1 < nl_1; ++l1)
            {
                const value_type facVal = vt->eval(l0, l1);
                out_messages[0][l0] = std::min(out_messages[0][l0], facVal + in_messages[1][l1]);
                out_messages[1][l1] = std::min(out_messages[1][l1], facVal + in_messages[0][l0]);
            }
        }
    }


    template<class VT>
    inline void l1_factor_to_variable_messages(
        const VT * vt,
        const label_type nl,
        const typename VT::value_type beta,
        const typename VT::value_type ** in_messages,
        typename VT::value_type ** out_messages
    ){
        using value_type = typename VT::value_type;
        if(beta>=0){
            std::copy(in_messages[1], in_messages[1]+nl, out_messages[1]);
            std::copy(in_messages[0], in_messages[0]+nl, out_messages[0]);
            // "forward pass"
            for(label_type l=1; l<nl; ++l){
                out_messages[0][l] = std::min(out_messages[0][l], out_messages[0][l-1] + beta);
                out_messages[1][l] = std::min(out_messages[1][l], out_messages[1][l-1] + beta);
            }
            // backward pass
            for(label_type l = nl-2; l>=0; --l){
                out_messages[0][l] = std::min(out_messages[0][l], out_messages[0][l+1] + beta);
                out_messages[1][l] = std::min(out_messages[1][l], out_messages[1][l+1] + beta);
                if(l == 0)
                    break;
            }
        }
        else{
            generic_second_order_factor_to_variable_messages(vt, nl, nl, in_messages, out_messages);
        }
    }


}






    template<class T>
    class XArrayTensor;

    template<class T, std::size_t NUM_LABELS>
    class StaticNumLabelTensor;

    template<class T>
    class TensorBase{
    public:
        using value_type = T;
        using label_type = std::size_t;
        using shape_type = arity_vector<label_type>;

        virtual ~TensorBase() {}
        virtual T operator[](const label_type * labels)const = 0;
        virtual T at(const label_type * labels_begin, const label_type * labels_end)const = 0;
        virtual std::size_t arity() const = 0;
        virtual std::size_t shape(const std::size_t)const = 0;
        virtual shape_type shape() const  = 0;
        virtual std::size_t sum_of_shape()const = 0;
        virtual void copy_corder(value_type * out)const = 0;
        virtual void add_values(value_type * out)const = 0;

        virtual std::unique_ptr<TensorBase<T>> clone()const = 0;

        virtual void factor_to_variable_messages(
            const value_type ** in_messages,
            value_type ** out_messages
        )const = 0;


        // virtual vbind()const = 0;
        // virtual ebind()const = 0;


        virtual std::unique_ptr<TensorBase<T>> bind(
            gsl::span<const std::size_t> positions,
            gsl::span<const label_type> labels
        )const = 0;

        // convert to a tensor with only binary labels but higher
        // arity.
        // undefined for tensors which are already binary
        virtual std::unique_ptr<TensorBase<T>> binarize()const = 0;
    };





    template<class T, class DERIVED>
    class TensorCrtpBase : public TensorBase<T>, public CrtpBase<DERIVED>
    {
    public:
        using value_type  = T;
        using shape_type = typename TensorBase<T>::shape_type;
        using label_type = typename TensorBase<T>::label_type;

        void factor_to_variable_messages(
            const value_type ** in_messages,
            value_type ** out_messages
        )const override{
            const auto shape = this->derived_cast().shape();
            const auto arity = shape.size();
            arity_vector<std::size_t> labels(arity, 0);

            for(auto ai=0; ai<arity; ++ai)
            {
                std::fill(out_messages[ai], out_messages[ai]+shape[ai], std::numeric_limits<value_type>::infinity());
            }

            detail::for_each_state(arity, shape, labels, [&](auto && lables){

                // sum of in fac(labels) + var-to-fac message at current labels
                auto e = this->derived_cast().operator[](labels.data());
                for(auto ai=0; ai<arity; ++ai)
                {
                    e += in_messages[ai][labels[ai]];
                }

                for(auto ai=0; ai<arity; ++ai)
                {
                    const auto label = labels[ai];
                    out_messages[ai][label] = std::min(
                        out_messages[ai][label], e - in_messages[ai][label]
                    );
                }
            });
        }

        value_type at(const label_type * labels_begin, const label_type * labels_end)const override{
            const auto arity = this->derived_cast().arity();
            if(std::distance(labels_begin, labels_end) != arity)
            {
                throw std::runtime_error("labels must match arity");
            }
            for(std::size_t i=0; i<arity; ++i){
                if(labels_begin[i] >= this->derived_cast().shape(i))
                {
                    throw std::runtime_error("label is out of bounds");
                }
            }
            return this->derived_cast()[labels_begin];
        }

        template<class ... ARGS, typename  = meta::all_integral<ARGS...>>
        value_type operator()(ARGS && ... args)const{

            std::array<label_type, sizeof ...(ARGS)> labels{
                label_type(std::forward<ARGS>(args)) ...
            };
            return this->derived_cast().operator[](labels.data());
        }

        shape_type shape() const override
        {
            const auto arity = this->derived_cast().arity();
            shape_type s(arity);
            for(std::size_t i=0; i<arity; ++i){
                s[i] = this->derived_cast().shape(i);
            }
            return s;
        }

        std::size_t sum_of_shape()const override
        {
            const auto shape = this->derived_cast().shape();
            return std::accumulate(shape.begin(), shape.end(), 0);
        }

        void add_values(value_type * out)const override
        {
            const auto shape = this->derived_cast().shape();
            const auto arity = shape.size();
            arity_vector<std::size_t> labels(arity, 0);
            auto i = 0;
            detail::for_each_state(arity, shape, labels, [&](auto && lables){
                 // do smth with state
                out[i] += this->derived_cast().operator[](labels.data());
                ++i;
            });
        }

        void copy_corder(value_type * out)const override
        {
            const auto shape = this->derived_cast().shape();
            const auto arity = shape.size();
            arity_vector<std::size_t> labels(arity, 0);
            auto i = 0;
            detail::for_each_state(arity, shape, labels, [&](auto && lables){
                 // do smth with state
                out[i] = this->derived_cast().operator[](labels.data());
                ++i;
            });
        }

        std::unique_ptr<TensorBase<T>> clone() const override{
            return std::make_unique<DERIVED>(this->derived_cast());
        }
        std::unique_ptr<TensorBase<T>> bind(
            gsl::span<const std::size_t> positions,
            gsl::span<const label_type> labels
        )const override {

            using binded_tensor_type = XArrayTensor<T>;
            using xshape_type = typename binded_tensor_type::xshape_type;

            const auto & self = this->derived_cast();
            const auto arity = self.arity();
            const auto sub_arity =  arity - positions.size();


            arity_vector<label_type> labels_buffer(arity);
            arity_vector<label_type> sub_labels_buffer(sub_arity);


            arity_vector<std::size_t> free_pos(sub_arity);

            // build shape of subtensor and
            // set the labels of the fixed variables
            xshape_type sub_shape(sub_arity);
            auto si=0;
            auto fixed_index = 0;
            for(auto ai=0; ai<arity; ++ai)
            {
                // todo this could be done with binary search
                const auto found = std::find(positions.begin(), positions.end(), ai) != positions.end();
                if(!found){
                    sub_shape[si] = self.shape(ai);
                    // free variable
                    free_pos[si] = ai;
                    ++si;
                } else{
                    // fixed variable
                    labels_buffer[ai] = labels[fixed_index];
                    ++fixed_index;
                }
            }

            auto tensor = std::make_unique<XArrayTensor<T>>(sub_shape);
            auto iter = tensor->xexpression().begin();

            // iterate over all states of the sub-tensor
            detail::for_each_state(sub_arity, sub_shape, sub_labels_buffer, [&](auto && lables){

                auto sub_ai = 0;
                for(auto i : free_pos)
                {
                    labels_buffer[i] = sub_labels_buffer[sub_ai];
                    ++sub_ai;
                }
                tensor->xexpression()[sub_labels_buffer] = self[labels_buffer.data()];
                // *iter = self[labels_buffer.data()];
                // ++iter;
            });

            return tensor;
        }

        std::unique_ptr<TensorBase<T>> binarize()const override{
            const auto & derived = this->derived_cast();
            const auto arity = derived.arity();
            const auto shape =  derived.shape();

            // how many binary variables are needed for
            // each of the variables
            arity_vector<std::size_t> num_binary_var(arity);

            // cumsum of num_binary_var st offset[0] == 0
            arity_vector<std::size_t> offset(arity, 0);

            // compute entries of `num_binary_var` and `offset`
            // and the actual arity of the binary tensor
            auto binary_arity = std::size_t(0);
            for(auto i=0; i<arity; ++i)
            {
                offset[i] = binary_arity;
                const auto n_bits =  static_cast<std::size_t>(std::ceil(std::log2(derived.shape(i))));
                num_binary_var[i] = n_bits;
                binary_arity += n_bits;
            }

            // allocate the result binary tensor
            auto binary_tensor = std::make_unique<StaticNumLabelTensor<value_type, 2>>(binary_arity);

            // fill the binary tensor

            // buffer to hold the labels for this tensor
            arity_vector<std::size_t> labels(arity);

            // buffer to hold the labels for the binary tensor
            binary_arity_vector<std::size_t> binary_labels(binary_arity);


            // iterate over all states of this tensor and translate
            // the labels st we have also the "binarized labels"
            detail::for_each_state(arity, shape, labels, [&](auto && lables){

                // the value / energy of the tensor itself
                // at the current label
                const auto val =  derived[labels.data()];

                // binarize the labels
                for(auto i=0; i<arity; ++i)
                {
                    auto binary_labels_begin_ptr = binary_labels.data() + offset[i];
                    detail::encode_binary(labels[i], binary_labels_begin_ptr, num_binary_var[i]);
                }

                // write the value of the tensor itself at the current label
                // to the output tensor at the appropriate binary labels
                binary_tensor->operator[](binary_labels.data()) =  val;
            });

            return std::move(binary_tensor);
        }

    };


    template<class T, std::size_t ARITY>
    class BinaryMultilinearTensor :  public TensorCrtpBase<T, BinaryMultilinearTensor<T, ARITY>>
    {
    public:
        using base_type = TensorBase<T>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        BinaryMultilinearTensor(const value_type beta)
        :
        m_beta(beta){
        }
        std::size_t sum_of_shape()const override{
            return ARITY * 2;
        }
        T operator[](const label_type * labels)const override{
            if(std::any_of(labels, labels+ARITY, [](auto l){return l == 0;}))
            {
                return 0.0;
            }
            return m_beta;
        }

        std::size_t arity()const override{
            return ARITY;
        }
        std::size_t shape(const std::size_t) const override{
            return 2;
        }
    private:
        value_type m_beta;
    };



    template<class T, std::size_t ARITY>
    class PottsNTensor :  public TensorCrtpBase<T, PottsNTensor<T, ARITY>>
    {
    public:
        using base_type = TensorBase<T>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        PottsNTensor(const std::size_t num_labels,const value_type beta)
        :   m_num_labels(num_labels),
            m_beta(beta){
        }
        std::size_t sum_of_shape()const override{
            return ARITY * m_num_labels;
        }
        T operator[](const label_type * labels)const override{
            const auto l0 = labels[0];
            for(auto i=1; i<ARITY;++i){
                if(labels[i] != l0){
                    return m_beta;
                }
            }
            return 0.0;
        }

        std::size_t arity()const override{
            return ARITY;
        }
        std::size_t shape(const std::size_t) const override{
            return m_num_labels;
        }

    private:
        std::size_t m_num_labels;
        value_type m_beta;
    };


    template<class T>
    class Potts2Tensor : public TensorCrtpBase<T, Potts2Tensor<T>>
    {
        //: public TensorBase<T>{
    public:
        using base_type = TensorCrtpBase<T, Potts2Tensor<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        Potts2Tensor(const std::size_t num_labels = 0,const value_type beta = value_type(0))
        :   m_num_labels(num_labels),
            m_beta(beta){
        }
        std::size_t sum_of_shape()const override{
            return 2 * m_num_labels;
        }
        virtual T operator[](const label_type * labels)const override{
            return labels[0] == labels[1] ? value_type(0) : m_beta;
        }
        virtual std::size_t arity()const override{
            return 2;
        }
        virtual std::size_t shape(const std::size_t) const override{
            return m_num_labels;
        }


        void factor_to_variable_messages(
            const value_type ** in_messages,
            value_type ** out_messages
        )const override{
            detail::potts2_factor_to_variable_messages(m_num_labels, m_beta, in_messages, out_messages);
        }
    private:
        std::size_t m_num_labels;
        value_type m_beta;
    };



    template<class T>
    class L1Tensor : public TensorCrtpBase<T, L1Tensor<T>>
    {
        //: public TensorBase<T>{
    public:
        using base_type = TensorCrtpBase<T, L1Tensor<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        L1Tensor(const std::size_t num_labels = 0,const value_type beta = value_type(0))
        :   m_num_labels(num_labels),
            m_beta(beta){
        }
        std::size_t sum_of_shape()const override{
            return 2 * m_num_labels;
        }
        T operator[](const label_type * labels)const override{
            return m_beta * std::abs(labels[0] - labels[1]);
        }
        std::size_t arity()const override{
            return 2;
        }
        std::size_t shape(const std::size_t) const override{
            return m_num_labels;
        }
        void factor_to_variable_messages(
            const value_type ** in_messages,
            value_type ** out_messages
        )const override{
            l1_factor_to_variable_messages(this, m_num_labels, m_beta, in_messages, out_messages);
        }
    private:
        std::size_t m_num_labels;
        value_type m_beta;
    };





    // template<class T>
    // class OneToOneTensor : public TensorCrtpBase<T, OneToOneTensor<T>>
    // {
    //     //: public TensorBase<T>{
    // public:
    //     using base_type = TensorCrtpBase<T, OneToOneTensor<T>>;
    //     using value_type = typename base_type::value_type;
    //     using label_type = typename base_type::label_type;
    //     using base_type::shape;

    //     OneToOneTensor(const std::size_t arity = 0, const std::size_t num_labels, const value_type beta = value_type(0))
    //     :   m_arity(arity),
    //         m_num_labels(num_labels),
    //         m_beta(beta),
    //         m_mutex(),
    //         m_used(arity)
    //     {
    //     }
    //     std::size_t sum_of_shape()const override{
    //         return m_arity * m_num_labels;
    //     }
    //     T operator[](const label_type * labels)const override{

    //         std::scoped_lock lock{m_mutex};
    //         std::fill(m_used.begin(), m_used.end(), 0);
    //         for(auto i=0; i<m_arity; ++i)
    //         {
    //             const auto l = labels[i];
    //             if(m_used[l] == 0)
    //             {
    //                 m_used[l] = 1;
    //             }
    //             else
    //             {
    //                 return m_beta;
    //             }
    //         }
    //         return 0.0;
    //     }
    //     std::size_t arity()const override{
    //         return m_arity;
    //     }
    //     std::size_t shape(const std::size_t) const override{
    //         return m_num_labels;
    //     }

    // private:
    //     std::size_t m_arity;
    //     std::size_t m_num_labels;
    //     value_type m_beta;
    //     std::mutex m_mutex;
    //     std::vector<uint8_t> m_used;
    // };


    template<class T>
    class DeltaUnary : public TensorCrtpBase<T, DeltaUnary<T>>
    {
        //: public TensorBase<T>{
    public:
        using base_type = TensorCrtpBase<T, DeltaUnary<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        DeltaUnary(const std::size_t num_labels=0, const  std::size_t label=0, value_type beta=0)
        :   m_num_labels(num_labels),
            m_label(label),
            m_beta(beta)
        {
        }
        constexpr std::size_t sum_of_shape()const override{
            return m_label;
        }
        T operator[](const label_type * labels)const override{
            return labels[0] == m_label ? value_type(0) : m_beta;
        }
        constexpr std::size_t arity()const override{
            return 1;
        }
        constexpr std::size_t shape(const std::size_t) const override{
            return m_label;
        }


    private:
        std::size_t m_num_labels;
        std::size_t m_label;
        value_type m_beta;
    };




    // optimized unary tensor for binary labels where 
    // we normalize the tenor st. f(0) = f*(0) - f*(1), f(1) = 0
    // => we only need to store a single value := f*(0) - f*(1)
    template<class T>
    class OptimizedBinaryUnary : public TensorCrtpBase<T, OptimizedBinaryUnary<T>>
    {
        //: public TensorBase<T>{
    public:
        using base_type = TensorCrtpBase<T, OptimizedBinaryUnary<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;
        using base_type::shape;

        OptimizedBinaryUnary(value_type val0 = value_type(0), value_type val1 = value_type(0))
        :   m_val0(val0-val1){
        }
        constexpr std::size_t sum_of_shape()const override{
            return 2;
        }
        T operator[](const label_type * labels)const override{
            return labels[0] == 0 ? m_val0 : value_type(0);
        }
        constexpr std::size_t arity()const override{
            return 1;
        }
        constexpr std::size_t shape(const std::size_t) const override{
            return 2;
        }

        void add_values(value_type * out)const override
        {
            out[0] += m_val0;
        }

        void copy_corder(value_type * out)const override
        {
            out[0] = m_val0;
            out[1] = value_type(0);
        }
    private:
        value_type m_val0;
    };





    template<class T>
    class UnaryTensor : public TensorCrtpBase<T, UnaryTensor<T>>
    {
    public:
        using base_type = TensorCrtpBase<T, UnaryTensor<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;

        using base_type::shape;

        template<class ITER>
        UnaryTensor(ITER values_begin, ITER values_end)
        :   m_values(values_begin, values_end)
        {
        }

        UnaryTensor(std::initializer_list<T> values)
        :   m_values(values.begin(), values.end())
        {
        }
        std::size_t sum_of_shape()const override{
            return m_values.size();
        }
        UnaryTensor(const label_type num_labels = 0, const value_type value = value_type(0))
        :   m_values(num_labels, value)
        {
        }

        T operator[](const label_type * labels)const override{
            return m_values[labels[0]];
        }
        T & operator[](label_type index){
            return m_values[index];
        }
        T operator[](label_type index)const{
            return m_values[index];
        }
        auto data() {
            return m_values.data();
        }
        auto data() const{
            return m_values.data();
        }
        std::size_t arity()const override{
            return 1;
        }
        std::size_t shape(const std::size_t) const override{
            return m_values.size();
        }

        void add_values(value_type * out)const override
        {
            for(auto i=0; i<m_values.size(); ++i)
            {
                out[i] += m_values[i];
            }
        }

        void copy_corder(value_type * out)const override
        {
            std::copy(m_values.begin(), m_values.end(), out);
        }
    private:
        std::vector<T> m_values;
    };






    template<class T, std::size_t NUM_LABELS>
    class StaticNumLabelTensor : public TensorCrtpBase<T, StaticNumLabelTensor<T, NUM_LABELS>>
    {
        //: public TensorBase<T>{
    public:
        using base_type = TensorCrtpBase<T, StaticNumLabelTensor<T, NUM_LABELS>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;

        using base_type::shape;

        StaticNumLabelTensor(const std::size_t arity)
        :   m_arity(arity),
            m_values(nullptr)
        {
            m_values = new value_type[detail::ipow(NUM_LABELS, arity)];
        }
        ~StaticNumLabelTensor(){
            delete[] m_values;
        }
        std::size_t sum_of_shape()const override{
            return m_arity * NUM_LABELS;
        }


        value_type operator[](const label_type * labels)const override{
            return m_values[get_offset(labels)];
        }

        value_type & operator[](const label_type * labels){
            return m_values[get_offset(labels)];
        }

        template<class ... ARGS, typename meta::all_integral<ARGS ...>>
        value_type & operator()(ARGS && ... args){

            std::array<label_type, sizeof ...(ARGS)> labels{
                label_type(std::forward<ARGS>(args)) ...
            };
            return this->derived_cast().operator[](labels.data());
        }


        std::size_t arity()const override{
            return m_arity;
        }
        constexpr std::size_t shape(const std::size_t) const override{
            return NUM_LABELS;
        }

    private:

        auto get_offset(const label_type * labels)const{
            auto stride = 1;
            auto offset = 0;
            for(auto i=m_arity;  i !=0; --i){
                offset += labels[i-1] * stride;
                stride *= NUM_LABELS;
            }
            return offset;
        }
        value_type * m_values;
        std::size_t m_arity;
    };


    template<class T>
    class XArrayTensor  : public TensorCrtpBase<T, XArrayTensor<T>>
    {
    public:
        using base_type = TensorCrtpBase<T, XArrayTensor<T>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;


        using xarray_type = xt::xarray<T>;
        using xshape_type = typename xarray_type::shape_type;

        using base_type::shape;


        XArrayTensor()
        : m_xarray(){
        }

        template<typename... Args>
        XArrayTensor(Args && ... args)
        : m_xarray(std::forward<Args>(args) ...)
        {

        }

        T operator[](const label_type * labels)const override{
            const auto label_span = gsl::span<const label_type>(labels, m_xarray.dimension());
            return m_xarray[label_span];
        }
        std::size_t arity() const override{
            return m_xarray.dimension();
        }
        std::size_t shape(const std::size_t d)const override{
            return m_xarray.shape(d);
        }

        void copy_corder(value_type * out)const override{
            std::copy(m_xarray.begin(), m_xarray.end(), out);
        }
        void add_values(value_type * out)const override{
            std::for_each(m_xarray.begin(), m_xarray.end(),[&](auto val){
                *out = val;
                ++out;
            });
        }

        auto & xexpression(){
            return m_xarray;
        }
        const auto & xexpression()const{
            return m_xarray;
        }

    private:
        xarray_type m_xarray;
    };


    template<class T, std::size_t ARITY>
    class XTensorTensor  : public TensorCrtpBase<T, XTensorTensor<T, ARITY>>
    {
    public:
        using base_type = TensorCrtpBase<T, XTensorTensor<T, ARITY>>;
        using value_type = typename base_type::value_type;
        using label_type = typename base_type::label_type;


        using xtensor_type = xt::xtensor<T, ARITY>;
        using xshape_type = typename xtensor_type::shape_type;

        using base_type::shape;


        XTensorTensor()
        : m_xtensor(){
        }

        template<typename... Args>
        XTensorTensor(Args && ... args)
        : m_xtensor(std::forward<Args>(args) ...)
        {

        }

        T operator[](const label_type * labels)const override{
            const auto label_span = gsl::span<const label_type>(labels, m_xtensor.dimension());
            return m_xtensor[label_span];
        }
        std::size_t arity() const override{
            return ARITY;
        }
        std::size_t shape(const std::size_t d)const override{
            return m_xtensor.shape(d);
        }

        void copy_corder(value_type * out)const override{
            std::copy(m_xtensor.begin(), m_xtensor.end(), out);
        }
        void add_values(value_type * out)const override{
            std::for_each(m_xtensor.begin(), m_xtensor.end(),[&](auto val){
                *out = val;
                ++out;
            });
        }

        auto & xexpression(){
            return m_xtensor;
        }
        const auto & xexpression()const{
            return m_xtensor;
        }

    private:
        xtensor_type m_xtensor;
    };

}