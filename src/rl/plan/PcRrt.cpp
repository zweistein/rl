#include <iostream>
#include <cmath>
#include <ctime>

#include <boost/make_shared.hpp>
#include <boost/random.hpp>

#include <Eigen/Eigenvalues>

#include <rl/sg/Shape.h>
#include <rl/sg/solid/Scene.h>

#include "SimpleModel.h"
#include "Viewer.h"
#include "GaussianSampler.h"

#include "PcRrt.h"

namespace rl
{
  namespace plan
  {
    ::std::string PcRrt::getName() const 
    {
      return "PCRRT";
    }

    PcRrt::PcRrt() :
      Rrt()
    {
      // this->gen = ::boost::make_shared<boost::random::mt19937>(std::time(0));
      std::cout << "fixed seed" << std::endl;
      this->gen = ::boost::make_shared<boost::random::mt19937>(42);
    }

    bool PcRrt::solve() 
    {
      // add the start configuation
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(*this->start));
      // set initial state based on one particle (covariance is zero)
      ::std::vector<::rl::math::Vector> v;
      v.push_back(*this->start);
      // this->tree[0][this->begin[0]].covariance = ::rl::math::Matrix::Zero(this->model->getDof(), this->model->getDof());
      this->tree[0][this->begin[0]].gState = ::boost::make_shared<GaussianState>(v);


      
      timer.start();
      timer.stop();
      
      while (timer.elapsed() < this->duration)
      {
        // vector for Voronoi vertex selection
        ::rl::math::Vector chosenSample(this->model->getDof());
        // choose new vertex using Voronoi bias
        this->choose(chosenSample);
        Neighbor n = this->nearest(this->tree[0], chosenSample);
        Vertex chosenVertex = n.first;
        
        // choose angle depending on covariance
        // ::rl::math::Vector dir = this->sampleDirection(chosenVertex);
        // float angle = ::std::atan2(dir[1], dir[0]) - ::std::atan2(0, 1);
        // float angle = this->sampleAngle(chosenVertex);
        // boost::random::mt19937 gen(std::time(0));
        boost::random::uniform_real_distribution<> distr(0, 2*M_PI);
        float angle = distr(*this->gen);

        ::rl::math::Matrix particles;
        // sampleParticles will return false if the particle set is not useful
        if (this->sampleParticles(chosenVertex, angle, this->nrParticles, particles))
        {
          // visualize particles
          this->drawParticles(particles);
          // fit a gaussian to the particles
          Gaussian gaussian(particles);

          Vertex newVertex = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(gaussian.mean));
          // this->tree[0][newVertex].covariance = gaussian.covariance;
          this->tree[0][newVertex].gState = ::boost::make_shared<GaussianState>(particles);


          this->addEdge(chosenVertex, newVertex, this->tree[0]);

          // try to connect to goal
          Neighbor nearest;
          nearest.first = newVertex;
          nearest.second = this->model->transformedDistance(*this->tree[0][newVertex].q, *this->goal);

          PossibleGoal possibleGoal;
          possibleGoal.neighbor = nearest;
          possibleGoal.q = this->tryConnect(this->tree[0], nearest, *this->goal);

          ::rl::math::Matrix goalParticles;
          if (NULL != possibleGoal.q && this->areEqual(*possibleGoal.q, *this->goal) && this->sampleGoalParticles(newVertex, *this->goal, this->nrParticles, goalParticles)) {
            this->drawParticles(goalParticles);
            Vertex connected = this->addVertex(this->tree[0], possibleGoal.q);
            this->addEdge(possibleGoal.neighbor.first, connected, this->tree[0]);
            this->end[0] = connected;
            return true;
          }
          // int foo;
          // std::cin >> foo;
        }
        
        timer.stop();
      }
      
      return false;
    }

    // based on number of empty clusters when clustering with different k
    bool PcRrt::isUnimodal(const ::rl::math::Matrix& particles)
    {
      int maxK = 10;
      ::std::vector<::rl::math::Real> sse(maxK);
      int emptyClusterCount = 0;
      for (int k = 1; k <= maxK; ++k)
      {
        ::std::vector<::std::vector<::rl::math::Vector> > clusters(k);
        this->kMeans(particles, k, clusters);

        sse[k-1] = 0.0;
        // for all clusters
        for (int i = 0; i < k; ++i)
        {
          if (clusters[i].empty())
          {
            emptyClusterCount++;
          }
        }
      }

      if (emptyClusterCount == maxK*(maxK+1)/2 - maxK) {
        return true;
      } else {
        return false;
      }

      // if (sseRatio > 15) return false;

      return true;
    }

    VectorPtr PcRrt::tryConnect(Tree& tree, const Neighbor& nearest, const ::rl::math::Vector& chosen)
    {
      ::rl::math::Real distance = nearest.second;
      ::rl::math::Real step = distance;
      

      bool reached = false;
      
      if (step <= this->delta)
      {
        reached = true;
      }
      else
      {
        step = this->delta;
      }
      
      VectorPtr last = ::boost::make_shared< ::rl::math::Vector >(this->model->getDof());
      
      this->model->interpolate(*tree[nearest.first].q, chosen, step / distance, *last);
      
      if (NULL != this->viewer)
      {
//        this->viewer->drawConfiguration(*last);
      }
      
      this->model->setPosition(*last);
      this->model->updateFrames();
      
      if (this->model->isColliding())
      {
        // return NULL;
        return VectorPtr();
      }
      
      ::rl::math::Vector next(this->model->getDof());
      
      while (!reached)
      {
        distance = this->model->distance(*last, chosen);
        step = distance;
        
        if (step <= this->delta)
        {
          reached = true;
        }
        else
        {
          step = this->delta;
        }
        
        this->model->interpolate(*last, chosen, step / distance, next);
        
        if (NULL != this->viewer)
        {
//          this->viewer->drawConfiguration(next);
        }
        
        this->model->setPosition(next);
        this->model->updateFrames();
        
        if (this->model->isColliding())
        {
          break;
        }
        
        *last = next;
      }
      
      return last;
    }

    bool PcRrt::sampleParticles(const Vertex& start, float angle, int nrParticles, ::rl::math::Matrix& particles)
    {
      boost::random::normal_distribution<> distr(angle, this->angleVariance);

      particles.resize(nrParticles, this->model->getDof());

      int fails = 0;
      int rowIdx = 0;

      std::vector<::std::string> shapes1, shapes2;
      std::string shape1 = "";
      std::string shape2 = "";

      while (rowIdx < nrParticles)
      {
        // Particle nextStep(*this->tree[0][start].q);
        Particle nextStep(this->model->getDof());
        this->tree[0][start].gState->sample(nextStep);

        ::rl::math::Real noisyAngle = distr(*this->gen);
        ::rl::math::Real stepX = ::std::cos(noisyAngle) * this->delta;
        ::rl::math::Real stepY = ::std::sin(noisyAngle) * this->delta;

        int steps = 0;
        do
        {
          nextStep[0] += stepX;
          nextStep[1] += stepY;

          this->model->setPosition(nextStep);
          this->model->updateFrames();
          steps++;
        }
        while (!this->model->isColliding());

        if (steps > 1)
        {
          std::string s1, s2;
          this->solidScene->lastCollidingShapes(s1, s2);
          if (shape1 == "")
          {
            shape1 = s1;
            shape2 = s2;
          }
          else if (s1 != shape1 || s2 != shape2)
          {
            // not every collision is between the same two shapes
            return false;
          }

          particles.row(rowIdx) = nextStep;
          rowIdx++;
        }
        else
        {
          fails++;
          if (fails > nrParticles / 4)
          {
            return false;
          }
        }
      }
      return true;
    }

    bool PcRrt::sampleGoalParticles(const Vertex& start, ::rl::math::Vector& goal, int nrParticles, ::rl::math::Matrix& particles)
    {
      particles.resize(nrParticles, this->model->getDof());

      ::rl::math::Vector startVec = *this->tree[0][start].q;
      ::rl::math::Vector dir = goal - startVec;
      ::rl::math::Real angle = ::std::atan2(dir[1], dir[0]) - ::std::atan2(0, 1);

      ::rl::math::Vector nextStep(startVec);

      int stepsToGo = 0;
      ::rl::math::Real stepX = ::std::cos(angle) * this->delta/4;
      ::rl::math::Real stepY = ::std::sin(angle) * this->delta/4;
      while (!this->areEqual(nextStep, goal))
      {
        nextStep[0] += stepX;
        nextStep[1] += stepY;
        stepsToGo++;
      }

      boost::random::normal_distribution<> angleDistr(angle, this->angleVariance);
      boost::random::normal_distribution<> stepDistr(0.0, this->stepVariance);

      int rowIdx = 0;

      while (rowIdx < nrParticles)
      {
        // ::rl::math::Vector nextStep(startVec);
        Particle nextStep(this->model->getDof());
        this->tree[0][start].gState->sample(nextStep);

        ::rl::math::Real noisyAngle = angleDistr(*this->gen);
        ::rl::math::Real stepX = ::std::cos(noisyAngle) * this->delta/4;
        ::rl::math::Real stepY = ::std::sin(noisyAngle) * this->delta/4;
        for (int step = 0; step < stepsToGo; step++)
        {
          nextStep[0] += stepX + stepDistr(*this->gen);
          nextStep[1] += stepY + stepDistr(*this->gen);
        }
        particles.row(rowIdx) = nextStep;

        rowIdx++;
      }

      return true;
    }
    
    ::rl::math::Vector PcRrt::sampleDirection(Vertex& vertex)
    {
      // calculate eigenvectors
      ::Eigen::EigenSolver<::rl::math::Matrix> eig(this->tree[0][vertex].gState->covariance());
      ::rl::math::Matrix evecs = eig.eigenvectors().real();
      // first principal component
      ::rl::math::Vector pc = evecs.col(0);
      // std::cout << "cov:" << std::endl << this->tree[0][vertex].covariance << std::endl << "first: " << std::endl << evecs.col(0) << std::endl << "second: " << std::endl << evecs.col(1) << std::endl;

      return pc.normalized();
      // boost::random::mt19937 gen(std::time(0));
      // boost::random::uniform_real_distribution<> distr(0, 2*M_PI);
      // return distr(*this->gen);
    }

    void PcRrt::drawParticles(::rl::math::Matrix& particles) 
    {
      // std::cout << "drawParticles" << std::endl;
      // std::cout << particles << std::endl;
      for (int rowIdx = 0; rowIdx < particles.rows(); ++rowIdx)
      {
        this->model->setPosition(particles.row(rowIdx));
        this->model->updateFrames();
        const ::rl::math::Transform& t = this->model->forwardPosition();
        ::rl::math::Vector3 opPos;
        opPos[0] = t.translation().x();
        opPos[1] = t.translation().y();
        opPos[2] = t.translation().z();
        // std::cout << opPos << std::endl;
        this->viewer->drawSphere(opPos, 0.02);
      }
    }

    void PcRrt::kMeans(const ::rl::math::Matrix& data, const int k, ::std::vector<::std::vector<::rl::math::Vector> >& clusters)
    {
      // just to be sure
      assert(data.rows() >= k && clusters.size() == k);

      // store cluster means in here
      ::std::vector<::rl::math::Vector> means(k);
      // init means with data points
      for (int i = 0; i < k; ++i)
      {
        // means[i] = data.row(i);
        means[i] = ::rl::math::Vector(this->model->getDof());
        this->choose(means[i]);
      }
      // store assigned cluster for each data point to notice changes
      ::std::vector<int> curClusters(data.rows(), -1);

      bool changed = true;
      while (changed)
      {
        for (auto& c : clusters)
        {
          c.clear();
        }

        changed = false;
        // assign data to clusters
        for (int dataIdx = 0; dataIdx < data.rows(); ++dataIdx)
        {
          int bestCluster = 0;
          ::rl::math::Real bestDist = this->model->transformedDistance(data.row(dataIdx), means[bestCluster]);
          for (int i = 1; i < k; ++i)
          {
            ::rl::math::Real dist = this->model->transformedDistance(data.row(dataIdx), means[i]);
            if (dist < bestDist)
            {
              dist = bestDist;
              bestCluster = i;
            }
          }
          clusters[bestCluster].push_back(data.row(dataIdx));
          if (curClusters[dataIdx] != bestCluster)
          {
            // this point has changed its cluster
            curClusters[dataIdx] = bestCluster;
            changed = true;
          }
        }

        // update means, keep the old mean to reduce risk of empty clusters
        for (int i = 0; i < k; ++i)
        {
          for (const auto& p : clusters[i])
          {
            means[i] += p;
          }
          means[i] /= clusters[i].size()+1;
        }
      }
    }
  }
}