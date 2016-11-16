//
// Copyright (c) 2009, Markus Rickert
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

#include <boost/make_shared.hpp>
#include <boost/random.hpp>
#include <boost/graph/random.hpp>

#include <rl/math/Vector.h>

#include "Est.h"
#include "SimpleModel.h"

#include <iostream>
#include <math.h>

namespace rl
{
  namespace plan
  {
    Est::Est() : Rrt() {}

    ::std::string 
    Est::getName() const 
    {
      return "EST";
    }

    bool Est::solve() 
    {
      ::std::cout << "Est solve!" << ::std::endl;
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared< ::rl::math::Vector >(*this->start));
      
      boost::random::mt19937 gen;
      boost::random::uniform_real_distribution<> distr(0, 2*M_PI);

      timer.start();
      timer.stop();

      Vertex chosenVertex = this->begin[0];
      ::rl::math::Real stepSize = 0.01;

      ::rl::math::Vector chosenSample(this->model->getDof());
      
      while (timer.elapsed() < this->duration)
      {
        int steps = 0;
        ::rl::math::Vector nextStep(*this->tree[0][chosenVertex].q);
        
        float angle = distr(gen);
        ::rl::math::Real stepX = ::std::cos(angle) * stepSize;
        ::rl::math::Real stepY = ::std::sin(angle) * stepSize; 
        
        while (!this->model->isColliding())
        {
          nextStep[0] += stepX;
          nextStep[1] += stepY;

          this->model->setPosition(nextStep);
          this->model->updateFrames();
          steps++;
        }

        // check if we actually moved through some free space
        if (steps > 1) 
        {
          Vertex collision_vertex = this->addVertex(this->tree[0], ::boost::make_shared< ::rl::math::Vector >(nextStep));
          this->addEdge(chosenVertex, collision_vertex, this->tree[0]);          
          
          // try to connect the new vertex to the goal
          Neighbor nearest;
          nearest.first = collision_vertex;
          nearest.second = this->model->transformedDistance(*this->tree[0][collision_vertex].q, *this->goal);

          Vertex connected = this->connect(this->tree[0], nearest, *this->goal);

          if (NULL != connected)
          {
            if (this->areEqual(*this->tree[0][connected].q, *this->goal)) 
            {
              this->end[0] = connected;
              return true;
            }
          }
        }

        // reset model to initial position, which should be collision-free
        this->model->setPosition(*this->tree[0][this->begin[0]].q);
        this->model->updateFrames();

        // choose new vertex randomly
        // chosenVertex = ::boost::random_vertex(this->tree[0], gen);

        // choose new vertex using Voronoi bias
        this->choose(chosenSample);
        chosenVertex = this->nearest(this->tree[0], chosenSample).first;

        timer.stop();
      }
      
      return false;
    }
  }
}
