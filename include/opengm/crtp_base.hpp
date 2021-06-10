#pragma once

namespace opengm {


    template<class derived>
    class CrtpBase{
    public:
        auto & derived_cast(){
            return *static_cast<derived*>(this);
        }
        const auto & derived_cast()const{
            return *static_cast<const derived*>(this);
        }
    };

}