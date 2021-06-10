#pragma once

#include <type_traits>

namespace opengm::meta{



    class null_type{};

    template<class T>
    struct is_null_type : public std::is_same<null_type, T>{};

    template<typename... Ts>
    using all_integral = typename std::enable_if<std::conjunction<std::is_integral<std::decay_t<Ts>>...>::value>::type;


    template<class A, class B>
    using if_not_null_type_t = std::conditional_t< !meta::is_null_type<A>::value,
        A,B
    >;

}