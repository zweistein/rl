//
// Copyright (c) 2016, Felix Wolff
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#ifndef _RL_PLAN_PCRRT_H_
#define _RL_PLAN_PCRRT_H_

#include <vector>
#include <memory>
#include <numeric>
#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

#include <Eigen/LU>

#include "Rrt.h"

namespace rl
{
  namespace plan
  {    
    /**
     * Particle Collision RRTs
     */
    class PcRrt : public Rrt {
    public:
      PcRrt();
      ::std::string getName() const;

      int nrParticles;
      
    protected:
      bool isUnimodal(const ::rl::math::Matrix& particles);

      typedef ::rl::math::Vector Particle;
      typedef ::std::vector<Particle> ParticleSet;

      class Gaussian {
      public:
        Gaussian(const ::rl::math::Matrix& particles)
        {
          this->init(particles);
        }

        Gaussian(const ::std::vector<::rl::math::Vector>& particles)
        {
          ::rl::math::Matrix p;
          p.resize(particles.size(), particles[0].size());
          for (int i = 0; i < particles.size(); ++i)
          {
            p.row(i) = particles[i];
          }
          this->init(p);
        }

        bool isUnimodal()
        {
          return true;
        }

        ::rl::math::Real likelihood(const ::std::vector<::rl::math::Vector>& data)
        {
          // std::cout << this->mean << std::endl << std::endl << this->covariance << std::endl;
          // check if covariance is of full rank
          ::rl::math::Real det = this->covariance.determinant();
          assert(det != 0.0);
          // if (det == 0) 
          // {
          //   // not full rank
          //   return -1;
          // }

          ::rl::math::Real lh = 0.0;
          for (const auto& d : data)
          {
            // ::rl::math::Real p1 = (d - this->mean).transpose() * this->covariance.inverse() * (d - this->mean);
            // std::cout << "p1 = " << p1 << std::endl;
            // ::rl::math::Real p = exp(-0.5 * p1) / sqrt(pow(2*M_PI, 2) * det);
            ::rl::math::Real p = exp(-0.5*(d - this->mean).transpose()*this->covariance.inverse()*(d - this->mean)) / sqrt((2*M_PI*this->covariance).determinant());
            // assert(false);
            // std::cout << d << std::endl << "p: " << p << std::endl;
            lh += log(p);
          }
          // avgProb /= data.size();

          // ::rl::math::Real lh = 0.0;
          // for (const auto& d : data)
          // {
          //   lh += -0.5*log(det) - 0.5*(d - this->mean).transpose()*this->covariance.inverse()*(d - this->mean) - 0.5*d.size()*log(2*M_PI);            
          // }

          return lh;
        }

        ::rl::math::Vector mean;
        ::rl::math::Matrix covariance;

      private:
        void init(const ::rl::math::Matrix& particles)
        {
          this->mean = particles.colwise().mean();
          // substract mean
          ::rl::math::Matrix centered = particles.rowwise() - this->mean.transpose();
          // calculate sample covariance
          this->covariance = particles.transpose() * particles / (particles.rows()-1);
        }
      };

      struct PossibleGoal
      {
        VectorPtr q;
        Neighbor neighbor;
      };

      virtual bool solve();
      // ::boost::shared_ptr<ParticleSet> sampleParticles(const Vertex& start, float angle, int nrParticles);
      bool sampleParticles(const Vertex& start, float angle, int nrParticles, ::rl::math::Matrix& particles);
      ::rl::math::Vector sampleDirection(Vertex& vertex);
      void drawParticles(::rl::math::Matrix& particles);
      virtual VectorPtr tryConnect(Tree& tree, const Neighbor& nearest, const ::rl::math::Vector& chosen);
      void kMeans(const ::rl::math::Matrix& data, const int k, ::std::vector<::std::vector<::rl::math::Vector> >& clusters);

    private:
      ::boost::shared_ptr<boost::random::mt19937> gen;
    };

  }
}

#endif // _RL_PLAN_PCRRT_H_
