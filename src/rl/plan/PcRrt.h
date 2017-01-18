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
#include <map>
#include <memory>
#include <numeric>
#include <cmath>
#include <tuple>

#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

#include <rl/sg/solid/Scene.h>
#include <rl/sg/Scene.h>

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
      ::rl::math::Real goalEpsilon;
      ::rl::sg::solid::Scene *solidScene;
      
    protected:
      virtual bool solve();

      void sampleDirection(::rl::math::Vector& rd);
      bool sampleConnectParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles);
      bool sampleGuardedParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles);
      bool sampleSlidingParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles);
      bool sampleGoalParticles(const Neighbor& nearest, ::rl::math::Vector& goal, int nrParticles, ::std::vector<Particle>& particles);

      bool isEqualCollisionState(::rl::sg::solid::Scene::CollisionMap& first, ::rl::sg::solid::Scene::CollisionMap& second);

      double projectOnSurface(const ::rl::math::Vector3& point, const ::rl::math::Vector3& pointOnSurface, const ::rl::math::Vector3& normal, ::rl::math::Vector3& out);
      bool moveConfigOnSurface(const ::rl::math::Vector& config, const ::rl::math::Vector3& pointOnSurface, const ::rl::math::Vector3& normal, ::rl::math::Vector& out);
      int expand(const VertexBundle& nearest, const ::rl::math::Transform& nearest2, const ::rl::math::Transform& chosen, const ::rl::math::Real& distance, VertexBundle& expanded);
      const rl::sg::solid::Scene::CollisionMap&  getAllCollidingShapes();
      void getPath(VectorList& path);
      bool getNormal(const Vertex& vertex, ::rl::math::Vector& normal);
      void drawParticles(const ::std::vector<Particle>& particles);
      void drawEigenvectors(Gaussian& gaussian, ::rl::math::Real scale);
      void drawSurfaceNormal(::rl::math::Vector3& startPoint, ::rl::math::Vector3& normal, ::rl::math::Real scale = 1.0);
      virtual VectorPtr tryConnect(Tree& tree, const Neighbor& nearest, const ::rl::math::Vector& chosen);
      //void kMeans(const ::rl::math::Matrix& data, const int k, ::std::vector<::std::vector<::rl::math::Vector> >& clusters);

      typedef ::std::vector<Particle> ParticleSet;


      struct PossibleGoal
      {
        VectorPtr q;
        Neighbor neighbor;
      };


    private:
      ::boost::shared_ptr<boost::random::mt19937> gen;
    };
  }
}

#endif // _RL_PLAN_PCRRT_H_
