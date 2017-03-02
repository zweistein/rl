#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

#include <boost/make_shared.hpp>
#include <boost/random.hpp>

#include <Eigen/Eigenvalues>

#include <rl/sg/Shape.h>
#include <rl/sg/Body.h>
#include <rl/sg/solid/Scene.h>

#include "NoisyModel.h"
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
      Rrt(),
      nrParticles(20),
      gamma(0.5)
    {
      // comment in for random seed
      int seed = std::time(0);
      this->gen = ::boost::make_shared<boost::random::mt19937>(seed);
      std::cout<<"seed: "<<seed<<std::endl;
      //this->gen = ::boost::make_shared<boost::random::mt19937>(44);
    }

    //Needed for Guarded moves
    void PcRrt::sampleDirection(::rl::math::Vector& rd)
    {
      boost::random::normal_distribution<> distr(0, 1);
      int dim = rd.rows();
      double rdsum = 0;
      for(int i=0; i<dim; i++)
      {
          rd[i] = distr(*this->gen);
          rdsum += rd[i]*rd[i];
      }
      rdsum = sqrt(rdsum);
      for(int i=0; i<dim; i++)
      {
          rd[i] /= rdsum;
      }
    }

    PcRrt::Neighbor PcRrt::nearest(const Tree& tree, const ::rl::math::Vector& chosen, double& directionSigma)
    {
      struct NeighborCompare
      {
        bool operator()(const Neighbor& lhs, const Neighbor& rhs) const
        {
          return lhs.second < rhs.second;
        }
      };

      std::set<Neighbor, NeighborCompare> neighbors;

      for (VertexIteratorPair i = ::boost::vertices(tree); i.first != i.second; ++i.first)
      {
        ::rl::math::Real d = this->model->transformedDistance(chosen, *tree[*i.first].q);
        neighbors.insert(Neighbor(*i.first, d));
      }

      Neighbor bestNeighbor(nullptr, ::std::numeric_limits<::rl::math::Real>::max());

      int tested = 0;
      for (auto n = neighbors.begin(); n != neighbors.end() && tested < 10; ++n)
      {
        const auto& vertex = (*n).first;
        const auto& mean = tree[vertex].gState->configMean();
        const auto& particles = tree[vertex].gState->getParticles();
        const ::rl::math::Real distance = this->model->distance(mean, chosen);

        std::vector<Particle> movedParticles;
        for (int pIdx = 0; pIdx < particles.size(); ++pIdx)
        {
          ::rl::math::Vector motionNoise(this->model->getDof());
          this->model->sampleMotionError(motionNoise);

          auto particle = particles[pIdx];
          auto particleOffset = particle.config - mean;
          auto shiftedChosen = chosen + particleOffset;

          ::rl::math::Vector dest(this->model->getDof());
          this->model->interpolateNoisy(particle.config, shiftedChosen, distance, motionNoise, dest);

          Particle p(dest);
          movedParticles.push_back(dest);
        }

        GaussianState g(movedParticles);
        ::rl::math::Real distanceMetric =  this->gamma*g.configGaussian().eigenvalues().sum();
        ::rl::math::Real uncertaintyMetric =  g.configGaussian().covariance.trace();
        ::rl::math::Real metric =  this->gamma * uncertaintyMetric + (1.0 - this->gamma) * distanceMetric;
        if (metric < bestNeighbor.second)
        {
          bestNeighbor.first = vertex;
          bestNeighbor.second = metric;
          rl::math::Vector dir = chosen-g.configMean();
          dir.normalize();
          //Mahalanobis distance squared
          //Eigen::LLT<Eigen::MatrixXd> lltOfG(g.configCovariance());

          //directionSigma = dir.transpose()*lltOfG.solve(dir);
        }

        tested++;
      }


      return bestNeighbor;
    }

    bool PcRrt::solve()
    {
      // add the start configuation
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(*this->start));
      // set initial state based on one particle (covariance is zero)

      ::std::vector<Particle> v;
      ::rl::math::Vector initialError(this->model->getDof());

      for(int i=0; i< this->nrParticles; i++)
      {
        ::rl::math::Vector sample = *this->start;
        do
        {
          this->model->sampleInitialError(initialError);

          for(int j=0; j<this->model->getDof(); j++)
          {
            sample[j]+=initialError[j];
          }
          this->model->setPosition(sample);
          this->model->updateFrames();
        }while(this->model->isColliding());

        Particle p(sample);
        v.push_back(p);

      }




      this->tree[0][this->begin[0]].gState = ::boost::make_shared<GaussianState>(v);
      this->model->setPosition(*this->start);
      this->model->updateFrames();
      this->tree[0][this->begin[0]].t = ::boost::make_shared<::rl::math::Transform>(this->model->forwardPosition());


      timer.start();
      timer.stop();

      this->model->setPosition(*this->goal);
      this->model->updateFrames();
      ::rl::math::Transform goalT = this->model->forwardPosition();

      while (timer.elapsed() < this->duration)
      {
        // vector for Voronoi vertex selection
        ::rl::math::Vector chosenSample(this->model->getDof());


        boost::random::uniform_int_distribution<> goalDistr(0, 100);
        if(goalDistr(*this->gen) < 10)
          chosenSample = *this->goal;

        // choose new vertex using Voronoi bias
        else
          this->choose(chosenSample);

        double directionSigma;
        Neighbor n = this->nearest(this->tree[0], chosenSample, directionSigma);

        Vertex chosenVertex = n.first;

        ::std::vector<Particle> particles;
        ::rl::math::Vector3 slidingNormal;
        slidingNormal.setZero();




        // sample[...]Particles will return false if the particle set is not useful
        bool sampleResult = false;


        if( this->tree[0][n.first].gState->isInCollision())
        {
          // randomly decide to do a slide or not
          boost::random::uniform_01<boost::random::mt19937> doSlideDistr(*this->gen);
          //double val = doSlideDistr()/(sqrt(directionSigma)*this->goalEpsilon);
          bool doSlide = doSlideDistr() < this->gamma;
          //std::cout<<"doslide "<<val<<std::endl;

          if (doSlide)
          {
            bool doGuardedSlide = doSlideDistr() < this->gamma;
            sampleResult = this->sampleSlidingParticles(doGuardedSlide, n, chosenSample, this->nrParticles, particles, slidingNormal);
          }
          else
            sampleResult = this->sampleConnectParticles(n, chosenSample, this->nrParticles, false, particles);
        }
        else
        {
          // randomly decide to do a slide or not
          boost::random::uniform_01<boost::random::mt19937> doGuardDistr(*this->gen);
          //double val = doGuardDistr()/(sqrt(directionSigma)*this->goalEpsilon);
          bool doGuardedMove = doGuardDistr() < this->gamma;
          //std::cout<<"doguard "<<val<<std::endl;

          if (doGuardedMove)
            sampleResult = this->sampleGuardedParticles(n, chosenSample, this->nrParticles, particles);
          else
            sampleResult = this->sampleConnectParticles(n, chosenSample, this->nrParticles, false, particles);
        }


        if (sampleResult)
        {
          // visualize particles
          //this->drawParticles(particles);

          ::boost::shared_ptr<GaussianState> gaussianState = ::boost::make_shared<GaussianState>(particles, slidingNormal);
          Gaussian g = gaussianState->configGaussian();
          //this->drawEigenvectors(g, 1.0);

//          if(particles.begin()->contacts.size()>0)
//            this->drawSurfaceNormal(particles.begin()->contacts.begin()->point,particles.begin()->contacts.begin()->normal_env,0.1);


          //std::cout<<g.eigenvalues().sum()<<std::endl;
          //if(g.eigenvalues().sum()<0.1)
          {
          // add a new vertex and edge
          VectorPtr mean = ::boost::make_shared<::rl::math::Vector>(g.mean);
          //this->viewer->drawConfiguration(g.mean);
          Vertex newVertex = this->addVertex(this->tree[0], mean);

          this->tree[0][newVertex].gState = gaussianState;
          this->model->setPosition(g.mean);
          this->model->updateFrames();
          this->tree[0][newVertex].t = ::boost::make_shared<::rl::math::Transform>(this->model->forwardPosition());

          this->addEdge(chosenVertex, newVertex, this->tree[0]);

          ::rl::math::Real maxError = 0;
          for(int i=0; i<particles.size(); i++)
          {
            this->model->setPosition(particles[i].config);
            this->model->updateFrames();
            maxError = std::max(maxError,::rl::math::transform::distance(goalT, this->model->forwardPosition(), 0.0));
          }

          //std::cout << "reached goal with error: " << maxError << " (max allowed: " << this->goalEpsilon << ")" << std::endl;

          if (maxError < this->goalEpsilon)
          {
            // visualize goal connect step
            this->drawParticles(particles);
            this->end[0] = newVertex;
            return true;
          }


          // try to connect to goal
          Neighbor nearest;
          nearest.first = newVertex;
          nearest.second = this->model->transformedDistance(*this->tree[0][newVertex].q, *this->goal);
          PossibleGoal possibleGoal;
          possibleGoal.neighbor = nearest;

          possibleGoal.q = this->tryConnect(this->tree[0], nearest, *this->goal);

          // check if connect step reached the goal
          if (NULL != possibleGoal.q && this->areEqual(*possibleGoal.q, *this->goal)) {
            // sample particles for the connect step
            ::std::vector<Particle> goalParticles;
            if (this->sampleConnectParticles(nearest, *this->goal, this->nrParticles, true, goalParticles))
            {
              ::boost::shared_ptr<GaussianState> goalState = ::boost::make_shared<GaussianState>(goalParticles);

              ::rl::math::Real maxError = 0;
              for(int i=0; i<goalParticles.size(); i++)
              {
                this->model->setPosition(goalParticles[i].config);
                this->model->updateFrames();
                maxError = std::max(maxError,::rl::math::transform::distance(goalT, this->model->forwardPosition(), 0.0));
              }

              std::cout << "reached goal with error: " << maxError << " (max allowed: " << this->goalEpsilon << ")" << std::endl;

              if (maxError < this->goalEpsilon)
              {
                // visualize goal connect step
                this->drawParticles(goalParticles);
                Gaussian g = goalState->configGaussian();
                // add goal connect step to tree
                Vertex connected = this->addVertex(this->tree[0], possibleGoal.q);
                this->tree[0][connected].gState = goalState;
                this->model->setPosition(*possibleGoal.q);
                this->model->updateFrames();
                this->tree[0][connected].t = ::boost::make_shared<::rl::math::Transform>(this->model->forwardPosition());

                this->addEdge(possibleGoal.neighbor.first, connected, this->tree[0]);
                this->end[0] = connected;
                return true;
              }
            }
          }
        }
        }

        timer.stop();
      }

      return false;
    }

    /**
      This is just the connect function of the RRT, but it does not insert a node into the tree.
    */
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

      this->model->interpolate(tree[nearest.first].gState->getParticles().begin()->config, chosen, step / distance, *last);
      this->model->setPosition(*last);
      this->model->updateFrames();

      ::rl::math::Vector next(this->model->getDof());

      int steps = 0;
      bool inInitialCollision = true;

      ::rl::sg::solid::Scene::CollisionMap initColls, allColls;

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
        this->model->setPosition(next);
        this->model->updateFrames();

        this->model->isColliding();

        if(steps == 0)
        {
          initColls = this->getAllCollidingShapes();
          if(initColls.size() == 0)
            inInitialCollision = false;
        }

        allColls = this->getAllCollidingShapes();
        if(inInitialCollision)
        {
          inInitialCollision = isEqualCollisionState(initColls, allColls);
        }

        //stuck in initial collision
        if(inInitialCollision && steps > 3)
          break;

        //unexpected collision
        if(!inInitialCollision && allColls.size()>0)
          break;


        *last = next;
        steps++;
      }

      return last;
    }

    /**
            Computes a path from start to goal and incorporates motion error if requested.
    */

    void PcRrt::getPath(VectorList& path, int count)
    {
      path.clear();
      Vertex i = this->end[0];
      ::rl::math::Vector q(this->model->getDof());;
      while (i != this->begin[0])
      {      
        ::boost::shared_ptr<GaussianState> gs = this->tree[0][i].gState;

        int idx = count % gs->getParticles().size();
        path.push_front(gs->getParticles()[idx].config);
        i = ::boost::source(*::boost::in_edges(i, this->tree[0]).first, this->tree[0]);
      }

      ::boost::shared_ptr<GaussianState> gs = this->tree[0][i].gState;
      int idx = count % gs->getParticles().size();
      path.push_front(gs->getParticles()[idx].config);

      if (this->useMotionError)
      {
        std::cout << "Using motion error" << std::endl;
        ::rl::math::Vector error;
        this->model->sampleMotionError(error);
        for (int i = 0; i < this->model->getDof(); ++i)
        {
          std::cout << "error joint " << i<< ": " << error[i] << std::endl;
        }

        ::std::vector<::rl::math::Vector> pathVec;
        ::std::vector<::rl::math::Vector> noisyPathVec;

        // convert list to vector for easy access
        for (const auto& p : path)
        {
          pathVec.push_back(p);
        }

        noisyPathVec.push_back(pathVec[0]);
        // apply error to path
        // -2 due to the fact, that the last step will not always move into collision
        // and needs to be treated differently
        for (int i = 0; i < pathVec.size()-2; ++i)
        {
          ::rl::math::Vector from = pathVec[i];
          ::rl::math::Vector to = pathVec[i+1];
          ::rl::math::Vector dir = to - from;

          ::rl::math::Vector noisyFrom = noisyPathVec[i];
          ::rl::math::Vector noisyTo(this->model->getDof());


          for(int i=0; i<this->model->getDof(); i++)
          {
            // apply noise
            noisyTo[i] = noisyFrom[i] + (1+error[i])*dir[i];
          }
          noisyPathVec.push_back(noisyTo);
        }

        // write noisy path
        path.clear();
        for (const auto& p : noisyPathVec)
        {
          path.push_back(p);
        }

        // write error to file
        ::rl::math::Real error_output = this->model->distance(noisyPathVec[noisyPathVec.size()-1], *this->goal);
        std::ofstream resultFile;
        resultFile.open(this->getName() + std::string("_results.txt"), std::ios::out | std::ios::app);
        resultFile << error_output << std::endl;

        //          // escape from current collision
        //          this->model->setPosition(noisyFrom);
        //          this->model->updateFrames();
        //          while (this->model->isColliding())
        //          {
        //            noisyFrom[0] += stepX;
        //            noisyFrom[1] += stepY;
        //            this->model->setPosition(noisyFrom);
        //            this->model->updateFrames();
        //          }
        //          // move into next collision
        //          while (!this->model->isColliding())
        //          {
        //            noisyFrom[0] += stepX;
        //            noisyFrom[1] += stepY;

        //            this->model->setPosition(noisyFrom);
        //            this->model->updateFrames();
        //          }

        //          noisyPathVec.push_back(noisyFrom);
        //        }

        //        // connect step
        //        // move calculated distance, or until collision
        //        ::rl::math::Vector from = pathVec[pathVec.size()-2];
        //        ::rl::math::Vector to = pathVec[pathVec.size()-1];
        //        ::rl::math::Vector dir = to - from;
        //        ::rl::math::Vector noisyFrom = noisyPathVec[pathVec.size()-2];

        //        // calculate noisy angle
        //        ::rl::math::Real noisyAngle = ::std::atan2(dir[1], dir[0]) - ::std::atan2(0, 1) + angleError;

        //        // pre-calculated distance to goal
        //        ::rl::math::Real distance = this->model->distance(from, to);

        //        ::rl::math::Real stepX = ::std::cos(noisyAngle) * this->delta;
        //        ::rl::math::Real stepY = ::std::sin(noisyAngle) * this->delta;

        //        ::rl::math::Real movedDistance = 0.0;
        //        this->model->setPosition(noisyFrom);
        //        this->model->updateFrames();
        //        // escape from current collision
        //        while (this->model->isColliding())
        //        {
        //          noisyFrom[0] += stepX;
        //          noisyFrom[1] += stepY;nt

        //          this->model->setPosition(noisyFrom);
        //          this->model->updateFrames();
        //        }

        //        // move into next collision
        //        while (!this->model->isColliding() && movedDistance < distance)
        //        {
        //          noisyFrom[0] += stepX;
        //          noisyFrom[1] += stepY;

        //          movedDistance += this->delta;

        //          this->model->setPosition(noisyFrom);
        //          this->model->updateFrames();
        //        }

        //        noisyPathVec.push_back(noisyFrom);

        //        // write noisy path
        //        path.clear();
        //        for (const auto& p : noisyPathVec)
        //        {
        //          path.push_back(p);
        //        }

        //        // write error to file
        //        ::rl::math::Real error = this->model->distance(noisyPathVec[noisyPathVec.size()-1], *this->goal);
        //        std::ofstream resultFile;
        //        resultFile.open(this->getName() + std::string("_results.txt"), std::ios::out | std::ios::app);
        //        resultFile << error << std::endl;
      }
    }

    void PcRrt::writeOutCurrentPath(std::string& path_string)
    {
      std::stringstream path_ss;
      Vertex i = this->end[0];
      std::list< Vertex > states;
      while (i != this->begin[0])
      {
        states.push_front(i);
        i = ::boost::source(*::boost::in_edges(i, this->tree[0]).first, this->tree[0]);
      }

      states.push_front(i);

      for( std::list< Vertex >::iterator it = states.begin(); it != states.end(); it++)
      {
        ::boost::shared_ptr<GaussianState> gs = this->tree[0][*it].gState;
        ::rl::math::Vector qm = gs->configMean();
        for(int i=0; i<qm.rows(); i++)
          path_ss<<qm(i)<<"\t";
        int inContact = 0;
        if(gs->isInCollision())
          inContact = 1;
        path_ss<<inContact<<"\t";

        if(gs->isSlidingMove())
        {
          path_ss<<gs->getSlidingNormal()(0)<<"\t"<<gs->getSlidingNormal()(1)<<"\t"<<gs->getSlidingNormal()(2)<<"\t";
        }
        else
        {
          path_ss<<"0\t0\t0\t";
        }

        ::rl::math::Transform ee_t = *(this->tree[0][*it].t);
        ::rl::math::Quaternion q(ee_t.linear());
        path_ss<<ee_t(0,3)<<"\t"<<ee_t(1,3)<<"\t"<<ee_t(2,3)<<"\t"<< q.x()<<"\t"<<q.y()<<"\t"<<q.z()<<"\t"<<q.w();
        path_ss<<std::endl;


      }
      path_string = path_ss.str();
      //std::cout<<path_string<<std::endl;
    }

    bool PcRrt::isEqualCollisionState(::rl::sg::solid::Scene::CollisionMap& first, ::rl::sg::solid::Scene::CollisionMap& second)
    {
      if (first.size() != second.size())
      {
        return false;
      }
      else
      {
        for(::rl::sg::solid::Scene::CollisionMap::iterator it = first.begin(); it != first.end(); it++)
        {
          if(second.find(it->first) == second.end())
          {
            return false;
          }
        }

        return true;
      }
    }

    /**
      Samples a set of particles for a move through free-space to a target configuration (chosen)
    */
    bool PcRrt::sampleConnectParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, bool goalConnect, ::std::vector<Particle>& particles)
    {
      particles.clear();

      int pIdx = 0;

      std::string shape1 = "";
      std::string shape2 = "";


      this->model->setPosition(this->tree[0][nearest.first].gState->configMean());
      this->model->updateFrames();
      ::rl::math::Vector3 initPoint = this->model->forwardPosition().translation();


      while (pIdx < nrParticles)
      {
        // sample a starting point from the gaussian
        ::rl::math::Vector init(this->model->getDof());

        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init, pIdx);

        ::rl::math::Vector nextStep = init;

        // sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->configMean();

        ::rl::math::Vector initialError = init - mean;
        ::rl::math::Vector target = chosen + initialError;

        int steps = 0;
        bool reached = false;
        bool collision = false;
        bool inInitialCollision = true;
        ::rl::math::Real step = this->delta;
        ::rl::math::Real distance = this->model->distance(init, target);

        rl::sg::solid::Scene::CollisionMap allColls, initColls;

        do
        {
          if (step >= distance)
          {
            reached = true;
            step = distance;
          }

          this->model->interpolateNoisy(init, target,  step / distance, motionNoise, nextStep);

          if(!this->model->isValid(nextStep))
            return false;

          this->model->setPosition(nextStep);
          this->model->updateFrames();

          collision = this->model->isColliding();

          if(steps == 0)
          {
            initColls = this->getAllCollidingShapes();
            if(initColls.size() == 0)
              inInitialCollision = false;
          }

          if(inInitialCollision)
          {
            allColls = this->getAllCollidingShapes();
            inInitialCollision = isEqualCollisionState(initColls, allColls);
          }

          if(inInitialCollision && steps > 3)
            return false;

          step += delta;
          steps++;
        }
        while (inInitialCollision || (!collision && !reached));

        if(goalConnect)
        {
          // valid particle, store it
          Particle p(nextStep);
          particles.push_back(p);

          pIdx++;
        }
        else
        {
          if (steps > 1 && reached && !collision)
          {
            Particle p(nextStep);
            particles.push_back(p);
            pIdx++;
            if (shape1 == "")
            {
              // first sample, init contact state shapes
              shape1 = "no_collision";
              shape2 = "no_collision";
            }
            else if (shape1 != "no_collision" || shape2 != "no_collision")
            {
              return false;
            }

          }
          else if (steps > 1 && collision)
          {
            rl::math::Vector3 normal;
            rl::sg::solid::Scene::CollisionMap collMap =  this->solidScene->getLastCollisions();

            //Body shape must be sensor
            if(!collMap.begin()->second.isSensor)
              return false;



            if(collMap.size()>1)
            {
              //std::cout<<"Connect move ends in a state with two collisions - this should not happen!"<<std::endl;
              return false;
            }


            std::string s1=collMap.begin()->first.first;
            std::string s2=collMap.begin()->first.second;

            if(!this->solidScene->getCollisionSurfaceNormal(initPoint,collMap.begin()->second.commonPoint,s1,s2,normal))
            {
              //this->viewer->drawLine(initPoint,collMap.begin()->second.commonPoint);
              std::cout<<"failed to compute normal"<<std::endl;
              return false;
            }


            if (shape1 == "")
            {
              // first sample, init contact state shapes
              shape1 = s1;
              shape2 = s2;
            }
            else if (s1 != shape1 || s2 != shape2)
            {
              // not every collision is between the same two shapes, -> not one contact state
              return false;
            }

            // valid particle, store it
            std::vector<Contact> contacts;
            rl::math::Vector3 contactPoint = collMap.begin()->second.commonPoint;
            contacts.push_back(Contact(contactPoint,normal,s1,s2));
            Particle p(nextStep, contacts);
            particles.push_back(p);

            pIdx++;
          }
          else
          {
            // we moved just one step until the next collision, or collided in last step
            return false;
          }
        }
      }

      return true;
    }

    /**
      Samples a set of particles for a move through free-space into a contact state.
    */
    bool PcRrt::sampleGuardedParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles)
    {
      particles.clear();

      int pIdx = 0;

      std::string shape1 = "";
      std::string shape2 = "";

      this->model->setPosition(this->tree[0][nearest.first].gState->configMean());
      this->model->updateFrames();
      ::rl::math::Vector3 initPoint = this->model->forwardPosition().translation();
//#define RANDOM_DIRECTION
#ifdef RANDOM_DIRECTION
      ::rl::math::Vector dir(this->model->getDof());
      sampleDirection(dir);
#endif

      while (pIdx < nrParticles)
      {
        rl::math::Vector init(this->model->getDof());
        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init, pIdx);

        // sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->configMean();

        ::rl::math::Vector initialError = init - mean;

#ifdef RANDOM_DIRECTION
        ::rl::math::Vector target = init + dir;
#else
        ::rl::math::Vector target = chosen + initialError;

#endif

        ::rl::math::Vector nextStep = init;

        // move until collision
        int steps = 0;

        bool collision;
        bool inInitialCollision = true;

        ::rl::sg::solid::Scene::CollisionMap allColls, initColls;

        ::rl::math::Real step = delta;
        ::rl::math::Real distance = this->model->distance(init, target);
        do
        {
          this->model->interpolateNoisy(init, target,  step / distance, motionNoise, nextStep);
          if(!this->model->isValid(nextStep))
            return false;

          this->model->setPosition(nextStep);
          this->model->updateFrames();

          collision = this->model->isColliding();

          if(steps == 0)
          {
            initColls = this->getAllCollidingShapes();
            if(initColls.size() == 0)
              inInitialCollision = false;
          }

          if(inInitialCollision)
          {
            allColls = this->getAllCollidingShapes();
            inInitialCollision = isEqualCollisionState(initColls, allColls);
          }

          if(inInitialCollision && steps > 3)
            return false;

          step += delta;
          steps++;
        }
        while (!collision || inInitialCollision);

        if (steps > 1)
        {

          rl::sg::solid::Scene::CollisionMap collMap =  this->solidScene->getLastCollisions();

          //collision must be perceivable by sensor
          if(!collMap.begin()->second.isSensor)
          {
        	std::cout << "Collision was not perceived by sensor\n";
            return false;
          }
          rl::math::Vector3 normal;
          std::string s1=collMap.begin()->first.first;
          std::string s2=collMap.begin()->first.second;

          if(!this->solidScene->getCollisionSurfaceNormal(initPoint,collMap.begin()->second.commonPoint,s1,s2,normal))
          {
            //this->viewer->drawLine(initPoint,collMap.begin()->second.commonPoint);
            std::cerr<<"failed to compute normal\n"; // for s1: "<< s1.c_str() <<" s2: "<< s2.c_str() <<std::endl;
            return false;
          } else { std::cout << "Norm computed \n"; }

          if(collMap.size()>1)
          {
            std::cerr<<"Guarded move ends in a state with two collisions - this should not happen! \n begin: "<<collMap.begin()->first.first << ","<< collMap.begin()->first.second << " end: " << collMap.end()->first.first<< "," << collMap.end()->first.second<<std::endl;
            return false;
          } else { std::cout << "collision with only one surface \n"; }


          if (shape1 == "")
          {
            // first collision, init contact state shapes
            shape1 = s1;
            shape2 = s2;
          }
          else if (s1 != shape1 || s2 != shape2)
          {
            // not every collision is between the same two shapes -> not one contact state
            return false;
          }

          // valid particle, store it
          std::vector<Contact> contacts;
          rl::math::Vector3 contactPoint = collMap.begin()->second.commonPoint;
          contacts.push_back(Contact(contactPoint,normal,s1,s2));
          Particle p(nextStep, contacts);
          particles.push_back(p);

          pIdx++;
        }
        else
        {
          // we moved just one step until the next collision, this is not useful
          return false;
        }
      }

      return true;
    }

    double PcRrt::projectOnSurface(const ::rl::math::Vector3& point, const ::rl::math::Vector3& pointOnSurface, const ::rl::math::Vector3& normal, ::rl::math::Vector3& out)
    {

      double dist = (point-pointOnSurface).dot(normal);
      out = point-dist*normal;
      return dist;
    }

    bool PcRrt::moveConfigOnSurface(const ::rl::math::Vector& config, const ::rl::math::Vector3& pointOnSurface, const ::rl::math::Vector3& normal, const std::pair<std::string, std::string>& collPair, ::rl::math::Vector& out)
    {

      //Step 1: move end-effector on the line given by point an normal

      double dist = std::numeric_limits<double>::infinity();
      out = config;

      this->model->setPosition(out);
      this->model->updateFrames();
      this->model->updateJacobian();
      this->model->updateJacobianInverse();

      ::rl::math::Transform eeT = this->model->forwardPosition();
      ::rl::math::Transform goalT = eeT;
      ::rl::math::Vector3 goaltransl;
      dist = projectOnSurface(eeT.translation(),pointOnSurface, normal, goaltransl);

      int steps = 0;
      while(fabs(dist) > 0.01 && steps++<30)
      {
        goalT.translation() = goaltransl;

        rl::math::Vector6 tdot;
        ::rl::math::transform::toDelta(eeT, goalT, tdot);

        rl::math::Vector qdot;
        this->model->inverseVelocity(tdot, qdot);

        if(qdot.norm() > this->delta)
        {
          qdot.normalize();
          qdot*=this->delta;
        }

        out += qdot;
        this->model->setPosition(out);
        this->model->updateFrames();
        this->model->updateJacobian();
        this->model->updateJacobianInverse();

        eeT = this->model->forwardPosition();
        dist = projectOnSurface(eeT.translation(),pointOnSurface, normal, goaltransl);
      }

      if(fabs(dist)>0.01)
        return false;

      //Step 2: Move end-effector towards surface until colliding
      steps = 0;

      this->model->isColliding();
      rl::sg::solid::Scene::CollisionMap colls = this->getAllCollidingShapes();
      bool collision = (colls.find(collPair) != colls.end());

      //Store the configuration before projecting onto collision surface
      ::rl::math::Vector preProjection = out;

      while(!collision && steps++ < 10)
      {
        rl::math::Vector6 tdot;
        tdot.setZero();
        tdot(0) = -normal(0); tdot(1) = -normal(1); tdot(2) = -normal(2);

        rl::math::Vector qdot;
        this->model->inverseVelocity(tdot, qdot);
        qdot.normalize();
        qdot *= this->delta;
        out += qdot;
        this->model->setPosition(out);
        this->model->updateFrames();
        this->model->updateJacobian();
        this->model->updateJacobianInverse();
        this->model->isColliding();
        colls = this->getAllCollidingShapes();
        collision = (colls.find(collPair) != colls.end());

      }



      //There is no contact with the sliding surface - return the configuration from step 1.
      if(!collision)
      {
        out = preProjection;
        return false;
      }

      return true;
    }

    /**
      Samples a set of particles for a sliding move along a surface.
    */
    bool PcRrt::sampleSlidingParticles(bool guardedMove, const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles, ::rl::math::Vector3& slidingNormal)
    {
      if (!this->tree[0][nearest.first].gState->isInCollision())
      {
        // we must start from within a collision to slide
        return false;
      }

      particles.clear();
      Particle oldParticle = *(this->tree[0][nearest.first].gState->getParticles().begin());

      //Sample a sliding surface
      int slidingIdx = 0;
      if(oldParticle.contacts.size() > 1)
      {
          boost::random::uniform_int_distribution<> chooseSurfaceDistr(0, oldParticle.contacts.size()-1);
          slidingIdx = chooseSurfaceDistr(*this->gen);
      }

      //The pose of a random sample
      this->model->setPosition(chosen);
      this->model->updateFrames();
      ::rl::math::Vector3 chosenForwardPos = this->model->forwardPosition().translation();


      Contact slidingContact = oldParticle.contacts[slidingIdx];
      std::pair<std::string, std::string> slidingPair;
      slidingPair.first  = slidingContact.shape_robot;
      slidingPair.second = slidingContact.shape_env;     
      slidingNormal = slidingContact.normal_env;

      // We move the EE frame to the contact point. This way the Jacobian
      // can be used to project the contact point back on the sliding surface.
      // Caution: whenever this function exits, model->resetTool must be called!
      this->model->setPosition(oldParticle.config);
      this->model->updateFrames();
      //ee in world frame
      ::rl::math::Transform ee_world = this->model->forwardPosition();
      //contact frame - orientation is same as ee
      ::rl::math::Transform cp_world = ee_world;
      //Translation ist shifted to contact point
      cp_world.translation() = slidingContact.point;
      //CP in EE frame
      ::rl::math::Transform cp_ee = ee_world.inverse()*cp_world;
      this->model->updateTool(cp_ee);

      //Now take the chosen  expansion direction and project it on the sliding surface
      ::rl::math::Vector3 chosen_proj;
      double distToSurface = projectOnSurface(chosenForwardPos, slidingContact.point, slidingNormal, chosen_proj);



      rl::sg::solid::Scene::CollisionMap allColls, finalCollisions;

      //The sliding velocity -- same for all particles
      ::rl::math::Vector6 tdot;

      int pIdx = 0;
      while (pIdx < nrParticles)
      {
        // sample a starting point
        ::rl::math::Vector init(this->model->getDof());
        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init, pIdx);

        //Sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        int steps = 0;

        //This is where the particle actually is
        ::rl::math::Vector nextStepReal = init;
        //This is where the particle thinks it is
        ::rl::math::Vector nextStepVirtual = this->tree[0][nearest.first].gState->configMean();

        bool reached = false;

        while (!reached)
        {
          this->model->setPosition(nextStepVirtual);
          this->model->updateFrames();
          this->model->updateJacobian();
          this->model->updateJacobianInverse();

          if(!guardedMove || steps == 0)
          {
            ee_world = this->model->forwardPosition();

            //This is the target pose of our slide
            ::rl::math::Transform goal_world = ee_world;
            goal_world.translation() = chosen_proj;

            //Compute a 6d Velocity twist for the task-space controller
            ::rl::math::transform::toDelta(ee_world, goal_world, tdot);

            //2D models do not move in z-direction (this is a convention)
            if(this->model->getDof() <= 2)
              tdot(2) = 0;
          }

          ::rl::math::Vector qdot(this->model->getDof());
          qdot.setZero();
          this->model->inverseVelocity(tdot, qdot);

          if(qdot.norm() < this->delta)
            reached = true;
          else
          {
            qdot.normalize();
            qdot *= this->delta;
          }


          ::rl::math::Vector newStep(this->model->getDof());
          //The robot thinks it is exactly following the line
          nextStepVirtual = nextStepVirtual+qdot;

          //In reality there is noise
          this->model->interpolateNoisy(nextStepReal,nextStepReal+qdot,1,motionNoise,newStep);
          //Project back on sliding surface
          if(!moveConfigOnSurface(newStep, slidingContact.point, slidingContact.normal_env, slidingPair, nextStepReal))
          {
            //std::cout<<"Could not project cfg on surface!!"<<std::endl;
            break;
          }

          //Check for collisions and if model is valid
          if(!this->model->isValid(nextStepReal))
          {
            this->model->resetTool();
            return false;
          }

          //this->viewer->drawConfiguration(nextStepReal);


          steps++;

          if(steps > 6.2/this->delta)
          {
            this->model->resetTool();
            return false;
          }

          this->model->setPosition(nextStepReal);
          this->model->updateFrames();
          this->model->updateJacobian();
          this->model->updateJacobianInverse();

          //check for singularity
          if(this->model->getDof() > 3 && this->model->getManipulabilityMeasure()  < 1.0e-3f)
          {
            this->model->resetTool();
            return false;
          }

          //Check for collision
          this->model->isColliding();
          allColls = this->getAllCollidingShapes();

          if (allColls.find(slidingPair) != allColls.end())
          {
            // check if we are in free-space or if we have another collision
            if (allColls.size() > 1)
            {
              // there is another collision, so the sliding ends here
              break;
            }
          }
          else
          {
            // we lost contact
            break;
          }
        }


//        // some magic number
//        if (steps < 3)
//        {
//          // we did not move very far, this is considered failure
//          this->model->resetTool();
//          return false;
//        }

        this->model->setPosition(nextStepReal);
        this->model->updateFrames();
        this->model->isColliding();
        allColls = this->getAllCollidingShapes();

        if (pIdx == 0)
        {
          finalCollisions =  this->getAllCollidingShapes();
        }
        else
        {
          if(!this->isEqualCollisionState(allColls, finalCollisions))
          {
            //inconsistent (over particles) collision state at end of slide
            this->model->resetTool();
            return false;
          }
        }

        // valid particle, store it
        if(allColls.size() == 0)
        {
          Particle p(nextStepReal);
          particles.push_back(p);
        }
        else
        {
          std::vector<Contact> contacts;
          for(::rl::sg::solid::Scene::CollisionMap::iterator it = allColls.begin(); it != allColls.end(); it++)
          {
            //Body shape must be sensor
            if(!it->second.isSensor)
            {
              this->model->resetTool();
              return false;

            }


            std::string c1 = it->first.first;
            std::string c2 = it->first.second;

            ::rl::math::Vector3 contactNormal;
            if(c1==slidingPair.first && c2==slidingPair.second)
            {
              contactNormal = slidingNormal;
            }
            else
            {

              rl::math::Vector3 source = slidingContact.point+slidingNormal*0.1;
              rl::math::Vector3 target = it->second.commonPoint;

              if(!this->solidScene->getCollisionSurfaceNormal(source, target, c1, c2 , contactNormal))
              {
                //this->viewer->drawLine(source,target);
                this->model->resetTool();
                std::cout<<"failed to compute normal"<<std::endl;
                return false;
              }
            }

            contacts.push_back(Contact(it->second.commonPoint,contactNormal,c1,c2));
          }
          Particle p(nextStepReal, contacts);
          particles.push_back(p);
        }
        pIdx++;
      }
      this->model->resetTool();
      return true;
    }

// DEPRECATED

//    bool PcRrt::getNormal(const Vertex& vertex, ::rl::math::Vector& normal)
//    {
//      ::rl::math::Vector evals = this->tree[0][vertex].gState->configGaussian().eigenvalues();
//      ::rl::math::Matrix evecs = this->tree[0][vertex].gState->configGaussian().eigenvectors();

//      int minIdx, maxIdx;
//      evals.minCoeff(&minIdx);
//      evals.maxCoeff(&maxIdx);

//      normal = evecs.col(minIdx).normalized();
//      // magic number assures that there is sufficient uncertainty in one direction
//      return evals[maxIdx] > 0.001;
//    }

    const rl::sg::solid::Scene::CollisionMap&  PcRrt::getAllCollidingShapes()
    {
      return this->solidScene->getLastCollisions();
    }

    void PcRrt::drawParticles(const ::std::vector<Particle>& particles)
    {
      for (int i = 0; i < particles.size(); ++i)
      {
        ::rl::math::Vector3 opPos;

        if(particles[i].contacts.size() > 0)
        {
          for (int j = 0; j < particles[i].contacts.size(); ++j)
          {
            opPos = particles[i].contacts[j].point;
            if(this->viewer)
              this->viewer->drawSphere(opPos, 0.02);
          }
        }
        else
        {
          this->model->setPosition(particles[i].config);
          this->model->updateFrames();
          const ::rl::math::Transform& t = this->model->forwardPosition();

          opPos[0] = t.translation().x();
          opPos[1] = t.translation().y();
          opPos[2] = t.translation().z();
          if(this->viewer)
            this->viewer->drawSphere(opPos, 0.02);
        }

      }
    }

    void PcRrt::drawEigenvectors(Gaussian& gaussian, ::rl::math::Real scale)
    {
      ::rl::math::Vector3 start, end1, end2, end3;
      ::rl::math::Matrix vectors = gaussian.eigenvectors();
      ::rl::math::Vector values = gaussian.eigenvalues();
      ::rl::math::Vector eigen1 = vectors.col(0) * sqrt(values(0)) * scale;
      ::rl::math::Vector eigen2 = vectors.col(1) * sqrt(values(1)) * scale;


      this->model->setPosition(gaussian.mean);
      this->model->updateFrames();
      const ::rl::math::Transform& t1 = this->model->forwardPosition();
      start(0) = t1.translation().x();
      start(1) = t1.translation().y();
      start(2) = t1.translation().z();



      if(this->model->getDof()==2)
      {
        this->model->setPosition(gaussian.mean+eigen1);
        this->model->updateFrames();
        const ::rl::math::Transform& t2 = this->model->forwardPosition();
        end1(0) = t2.translation().x();
        end1(1) = t2.translation().y();
        end1(2) = t2.translation().z();

        this->model->setPosition(gaussian.mean+eigen2);
        this->model->updateFrames();
        const ::rl::math::Transform& t3 = this->model->forwardPosition();
        end2(0) = t3.translation().x();
        end2(1) = t3.translation().y();
        end2(2) = t3.translation().z();
      }

      if(this->model->getDof()==3)
      {
        ::rl::math::Vector eigen3 = vectors.col(2) * sqrt(values(2)) * scale;
        end1 = start+eigen1;
        end2 = start+eigen2;
        end3 = start+eigen3;
        this->viewer->drawLine(start, end3);
      }

      this->viewer->drawLine(start, end1);
      this->viewer->drawLine(start, end2);

    }

    void PcRrt::drawSurfaceNormal(::rl::math::Vector3& startPoint, ::rl::math::Vector3& normal, ::rl::math::Real scale)
    {
      this->viewer->drawLine(startPoint, startPoint+normal*scale);
    }

    //    void PcRrt::kMeans(const ::rl::math::Matrix& data, const int k, ::std::vector<::std::vector<::rl::math::Vector> >& clusters)
    //    {
    //      // just to be sure
    //      assert(data.rows() >= k && clusters.size() == k);

    //      // store cluster means in here
    //      ::std::vector<::rl::math::Vector> means(k);
    //      // init means with data points
    //      for (int i = 0; i < k; ++i)
    //      {
    //        // means[i] = data.row(i);
    //        means[i] = ::rl::math::Vector(this->model->getDof());
    //        this->choose(means[i]);
    //      }
    //      // store assigned cluster for each data point to notice changes
    //      ::std::vector<int> curClusters(data.rows(), -1);

    //      bool changed = true;
    //      while (changed)
    //      {
    //        for (auto& c : clusters)
    //        {
    //          c.clear();
    //        }

    //        changed = false;
    //        // assign data to clusters
    //        for (int dataIdx = 0; dataIdx < data.rows(); ++dataIdx)
    //        {
    //          int bestCluster = 0;
    //          ::rl::math::Real bestDist = this->model->transformedDistance(data.row(dataIdx), means[bestCluster]);
    //          for (int i = 1; i < k; ++i)
    //          {
    //            ::rl::math::Real dist = this->model->transformedDistance(data.row(dataIdx), means[i]);
    //            if (dist < bestDist)
    //            {
    //              dist = bestDist;
    //              bestCluster = i;
    //            }
    //          }
    //          clusters[bestCluster].push_back(data.row(dataIdx));
    //          if (curClusters[dataIdx] != bestCluster)
    //          {
    //            // this point has changed its cluster
    //            curClusters[dataIdx] = bestCluster;
    //            changed = true;
    //          }
    //        }

    //        // update means, keep the old mean to reduce risk of empty clusters
    //        for (int i = 0; i < k; ++i)
    //        {
    //          for (const auto& p : clusters[i])
    //          {
    //            means[i] += p;
    //          }
    //          means[i] /= clusters[i].size()+1;
    //        }
    //      }
    //    }
  }
}
