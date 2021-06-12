#include <doctest.h>
#include <vector>
#include <sstream>
#include <random>



#include "opengm/minimizer/brute_force_naive.hpp"



// the standard forbids writing in the std namespace but it works on all compilers
namespace std
{
template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& in) {
    stream << "[";
    for(size_t i = 0; i < in.size(); ++i)
        if(i < in.size() - 1)
            stream << in[i] << ", ";
        else
            stream << in[i];
    stream << "]";
    return stream;
}
}




namespace opengm{

template<class GM>
auto solve_brute_force(GM && gm){

    using gm_type = std::decay_t<GM>;
    using minimizer_type = opengm::BruteForceNaive<gm_type>;
    using settings_type = typename minimizer_type::settings_type;
    settings_type settings;
    minimizer_type minimizer(gm, settings);
    minimizer.minimize();
    return std::make_pair(minimizer.best_labels(),  gm.evaluate(minimizer.best_labels()));
}





namespace condition
{
    class Optimal{
    public:
        template<class MINIMIZER>
        void operator()(MINIMIZER & minimizer)const{
            const auto & gm = minimizer.gm();

            const auto best_labels = minimizer.best_labels();
            const auto best_energy = minimizer.best_energy();

            const auto [opt_labels,opt_energy] = solve_brute_force(gm);
            if(best_labels != opt_labels)
            {
                std::stringstream ss;
                ss <<"minimizer did not yield optimal results\n";
                ss<<"BEST    "<<best_labels<<"\n";
                ss<<"OPTIMAL "<<opt_labels<<"\n";
                CHECK_MESSAGE(best_energy <= opt_energy, ss.str().c_str());
            }
        };
    };

    class NoCondition{
    public:
        template<class MINIMIZER>
        void operator()(MINIMIZER & minimizer)const{
        }
    };

    template<std::size_t K>
    class KOptimal;

    template<>
    class KOptimal<1>
    {
    public:
        template<class MINIMIZER>
        void operator()(MINIMIZER & minimizer)const{
            const auto & gm = minimizer.gm();
            const auto best_labels = minimizer.best_labels();
            const auto best_energy = minimizer.best_energy();

            auto labels = best_labels;
            for(auto vi=0; vi<gm.num_variables(); ++vi)
            {
                const auto nl = gm.num_labels(vi);
                for(auto l=0; l<nl; ++l)
                {
                    if(l != best_labels[vi])
                    {
                        labels[vi] = l;
                        auto e = gm.evaluate(labels);
                        CHECK_MESSAGE(e >= best_energy, "labels are not 1-Optimal");

                    }
                }
                // restore best
                labels[vi] = best_labels[vi];
            }
        };
    };


    template<>
    class KOptimal<2>
    {
    public:
        template<class MINIMIZER>
        void operator()(MINIMIZER & minimizer)const{
            const auto & gm = minimizer.gm();
            const auto best_labels = minimizer.best_labels();
            const auto best_energy = minimizer.best_energy();

            auto labels = best_labels;
            for(auto vi0=0; vi0<gm.num_variables()-1; ++vi0)
            {
                for(auto vi1=vi0+1; vi1<gm.num_variables(); ++vi1)
                {
                    const auto nl0 = gm.num_labels(vi0);
                    const auto nl1 = gm.num_labels(vi1);

                    for(auto l0=0; l0<nl0; ++l0)
                    {
                        for(auto l1=0; l1<nl1; ++l1)
                        {
                            if(l0 != best_labels[vi0] && l1 != best_labels[vi1])
                            {
                                labels[vi0] = l0;
                                labels[vi1] = l1;
                                auto e = gm.evaluate(labels);
                                CHECK_MESSAGE(e >= best_energy, "labels are not 2-Optimal");
                            }
                        }
                    }
                    // restore best
                    labels[vi0] = best_labels[vi0];
                    labels[vi1] = best_labels[vi1];
                }
            }
        };
    };

    class FactorLocalOptimal
    {
    public:
        template<class MINIMIZER>
        void operator()(MINIMIZER & minimizer)const{
            const auto & gm = minimizer.gm();
            const auto best_labels = minimizer.best_labels();
            const auto best_energy = minimizer.best_energy();

            CHECK(best_energy == doctest::Approx(gm.evaluate(best_labels)));

            auto labels = best_labels;
            auto fi = 0;
            for(auto && factor : gm)
            {
                auto && variables = factor.variables();
                gm.space().for_each_state(variables.begin(),variables.end(), labels, [&](auto && labels){
                    bool any_diff = false;
                    for(auto var : variables)
                    {
                        if(labels[var] != best_labels[var])
                        {
                            any_diff = true;
                            break;
                        }
                    }
                    if(any_diff)
                    {
                        std::stringstream ss;
                        ss<<"factor "<<fi<<" arity "<<factor.arity()<<"\n";
                        INFO(ss.str());
                        auto e = gm.evaluate(labels);
                        INFO("SOLVER ",best_labels);
                        INFO("FOUND  ",labels);
                        CHECK_MESSAGE(e >= best_energy, "labels are not FactorLocalOptional");
                    }
                });
                // restore best state
                for(auto var : variables)
                {
                    labels[var] = best_labels[var];
                }
                ++fi;
            }
        };
    };


}


template<class GM_GEN, class F, class CONDITION>
void testing(GM_GEN && gen, F && f, CONDITION && condition, std::size_t n_runs){

    SUBCASE(gen.name().c_str())
    {
        INFO(("gm-generator: "+gen.name()));
        for(size_t i=0; i<n_runs; ++i)
        {
            auto seed = gen.seed();
            auto gm = gen();
            auto sed = gen.seed();
            auto minimizer =  f(gm);

            INFO((std::string("minimizer: ")+minimizer.name()+ std::string(" seed ") + std::to_string(seed)));
            minimizer.minimize();
            auto best_labels = minimizer.best_labels();
            auto best_energy = minimizer.best_energy();

            condition(minimizer);
        }
    }
}

}