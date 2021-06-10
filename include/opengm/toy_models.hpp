#pragma once

#include <random>
#include <algorithm>
#include <string>

#include <xtensor/xrandom.hpp>

#include "opengm/graphical_model.hpp"
#include "opengm/tensors.hpp"
#include "opengm/space.hpp"


#include <iostream>

namespace opengm{


namespace detail{

    class RandomPottsGridImpl{

    public:
        RandomPottsGridImpl(
            std::size_t nx,
            std::size_t ny,
            std::size_t n_labels,
            bool make_tree,
            std::size_t seed = 42
        )
        :
        m_nx(nx),
        m_ny(ny),
        m_n_labels(n_labels),
        m_make_tree(make_tree),
        m_seed(seed){

        }


        auto seed(){
            return m_seed;
        }
        auto operator()(){


            using value_type = float;
            using value_ilist = std::initializer_list<value_type>;
            using label_type = std::size_t;
            using space_type = opengm::UniformSpace<label_type>;
            using GmType = opengm::GraphicalModel<space_type, value_type>;

            auto n_variables = m_nx*m_ny;
            auto get_vi = [&](auto x,auto y){
                return x*m_ny + y;
            };
            GmType gm(n_variables, m_n_labels);

            std::uniform_real_distribution<value_type> distribution(0.0f, 1.0f);
            std::default_random_engine generator(m_seed);
            ++m_seed;

            std::vector<value_type> data(m_n_labels);

            for(auto vi=0; vi<n_variables; ++vi){
                std::generate(data.begin(), data.end(), [&]() { return distribution(generator);});
                auto tensor = std::make_unique<opengm::UnaryTensor<value_type>>(data.begin(), data.end());
                gm.add_unary_factor(std::move(tensor), vi);
            }
            for(auto x=0; x<m_nx; ++x)
            {
                for(auto y=0; y<m_ny; ++y)
                {
                    auto vi = get_vi(x, y);
                    if( x+1 < m_nx)
                    {
                        if(!m_make_tree || y == 0){
                            auto tensor = std::make_unique<opengm::Potts2Tensor<value_type>>(m_n_labels, 0.1 * distribution(generator));
                            gm.add_factor(std::move(tensor), {vi,  get_vi(x+1, y)});
                        }

                    }
                    if(y+1 < m_ny)
                    {
                        auto tensor = std::make_unique<opengm::Potts2Tensor<value_type>>(m_n_labels, 0.1 * distribution(generator));
                        gm.add_factor(std::move(tensor), {vi,  get_vi(x, y+1)});

                    }
                }
            }
            //std::cout<<"seed "<<m_seed<<"\n";
            return gm;
        }



        std::string name()const{return "RandomPottsGridImpl";}
    private:
        std::size_t m_nx;
        std::size_t m_ny;
        std::size_t m_n_labels;
        bool m_make_tree;
        std::size_t m_seed;
    };



}

    struct RandomPottsGrid : public detail::RandomPottsGridImpl{
        RandomPottsGrid(
            std::size_t nx,
            std::size_t ny,
            std::size_t n_labels,
            std::size_t seed = 42
        ): detail::RandomPottsGridImpl(nx,ny,n_labels,false, seed){

        }
        std::string name()const{return "RandomPottsGrid";}
    };

    struct RandomPottsChain : public detail::RandomPottsGridImpl{
        RandomPottsChain(
            std::size_t n_variables,
            std::size_t n_labels,
            std::size_t seed = 42
        ): detail::RandomPottsGridImpl(n_variables,1,n_labels,false, seed){

        }
        std::string name()const{return "RandomPottsChain";}
    };

    struct RandomPottsGridFan : public detail::RandomPottsGridImpl{
        RandomPottsGridFan(
            std::size_t nx,
            std::size_t ny,
            std::size_t n_labels,
            std::size_t seed = 42
        ): detail::RandomPottsGridImpl(nx,ny,n_labels,true, seed){

        }
        std::string name()const{return "RandomPottsGridFan";}
    };


    template<class T = float>
    class RandomModel{

    public:

        using value_type = T;
        using value_ilist = std::initializer_list<value_type>;
        using label_type = std::size_t;
        using space_type = opengm::ExplicitSpace<label_type>;
        using GmType = opengm::GraphicalModel<space_type, value_type>;


        RandomModel(
            std::size_t n_var,
            std::size_t n_factors,
            std::size_t min_num_labels,
            std::size_t max_num_labels,
            std::size_t min_arity,
            std::size_t max_arity,
            std::size_t seed = 42
        )
        :   m_n_var(n_var),
            m_n_factors(n_factors),
            m_min_num_labels(min_num_labels),
            m_max_num_labels(max_num_labels),
            m_min_arity(min_arity),
            m_max_arity(max_arity),
            m_seed(seed)
        {

        }

        auto seed(){
            return m_seed;
        }
        auto operator()(){

            std::uniform_int_distribution<std::size_t> arity_distribution(m_min_arity, m_max_arity);
            std::uniform_int_distribution<std::size_t> num_labels_distribution(m_min_num_labels, m_max_num_labels);
            std::uniform_real_distribution<value_type> value_distribution(0.0f, 1.0f);
            std::default_random_engine generator(m_seed);
            ++m_seed;


            // space
            std::vector<label_type> spacevec(m_n_var);
            std::generate(spacevec.begin(), spacevec.end(), [&]() { return num_labels_distribution(generator);});

            // gm
            GmType gm(spacevec.begin(), spacevec.end());

            // all indices
            std::vector<std::size_t> vis(m_n_var);
            std::iota(vis.begin(), vis.end(), 0);

            // add factors
            for(auto i=0; i<m_n_factors; ++i)
            {
                // random arity of this factor
                const auto arity = arity_distribution(generator);

                // build random indices
                std::shuffle(vis.begin(), vis.end(), generator);

                // tensor shape
                std::vector<label_type> factor_shape(arity);
                for(auto d=0; d<arity; ++d)
                {
                    factor_shape[d] = gm.num_labels(vis[d]);
                }

                // tensor with random values
                auto tensor =  std::make_unique<XArrayTensor<value_type>>(
                    xt::random::rand(factor_shape, -1.0, 1.0)
                );

                gm.add_factor(std::move(tensor), vis.begin(), vis.begin() + arity);
            }



            return gm;
        }

        std::size_t m_n_var;
        std::size_t m_n_factors;
        std::size_t m_min_num_labels;
        std::size_t m_max_num_labels;
        std::size_t m_min_arity;
        std::size_t m_max_arity;
        std::size_t m_seed;
    };


}