//
// Copyright (c) 2016, Arne Sieverling
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

#include "NoisyModel.h"

namespace rl
{
  namespace plan
  {
    NoisyModel::NoisyModel() :
      DistanceModel()
    {
    }
    
    NoisyModel::~NoisyModel()
    {
    }

    void
    NoisyModel::interpolateNoisy(const ::rl::math::Vector& q1, const ::rl::math::Vector& q2, const ::rl::math::Real& alpha, const ::rl::math::Vector& noise, ::rl::math::Vector& q) const
    {

      // interpolate(q1,q2,alpha,q);
      // q=q+alpha*(q2-q1).cwiseProduct(noise);

      // noise[0] -> angle error
      // noise[1] -> step error

      ::rl::math::Vector direction = q2 - q1;
      ::rl::math::Real distance = this->distance(q1, q2) * alpha;
      distance += distance * noise[1];
      
      ::rl::math::Real angle = ::std::atan2(direction[1], direction[0]) - ::std::atan2(0, 1) + noise[0];
      ::rl::math::Real stepX = ::std::cos(angle) * distance;
      ::rl::math::Real stepY = ::std::sin(angle) * distance;

      q[0] = q1[0] + stepX;
      q[1] = q1[1] + stepY;
    }

    void
    NoisyModel::sampleInitialError(::rl::math::Vector &error)
    {
      if (NULL == this->motionErrorGen)
      {
        this->motionErrorGen = ::boost::make_shared<::boost::random::mt19937>(42);
      }

      if(this->motionError->rows()!=this->getDof())
      {
        std::cout << "warning: did not set initial error - will use default value 0.05" << std::endl;
        this->initialError->setOnes(this->getDof());
        (*this->motionError)*=0.05;
      }

      for(int i=0; i<this->getDof(); i++)
      {
        // sample an initial error
        ::boost::random::normal_distribution<> errorDistr(0, (*this->initialError)(i));
        error[i] = errorDistr(*this->motionErrorGen);
      }


    }

    void
    NoisyModel::sampleMotionError(::rl::math::Vector &error)
    {
      if (NULL == this->motionErrorGen)
      {
        this->motionErrorGen = ::boost::make_shared<::boost::random::mt19937>(42);
      }

      error.resize(this->getDof());

      if(this->motionError->rows()!=this->getDof())
      {
        std::cout << "warning: did not set motion error - will use default value 0.1" << std::endl;
        this->motionError->setOnes(this->getDof());
        (*this->motionError)*=0.1;
      }

      for(int i=0; i<this->getDof(); i++)
      {
        // sample a step error
        ::boost::random::normal_distribution<> stepDistr(0, (*this->motionError)(i));
        error[i] = stepDistr(*this->motionErrorGen);
        //std::cout << "error joint "<<i<<": " << error[i] << std::endl;
      }

    }
  }
}
