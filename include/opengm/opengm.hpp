#pragma once
#ifndef OPENGM_OPENGM_HPP
#define OPENGM_OPENGM_HPP

#include <cstdint>
#include <iostream>

namespace opengm {
    class MyClass
    {
    public:
        MyClass(const uint64_t size)
        : m_size(size)
        {

        }
        void hello_world()
        {
            std::cout<<"Hello World!\n";
        }
    private:
        uint64_t m_size;
    };

} // end namespace opengm


#endif // OPENGM_OPENGM_HPP