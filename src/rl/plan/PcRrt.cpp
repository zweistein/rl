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
      Rrt()
    {
      // comment in for random seed
      // this->gen = ::boost::make_shared<boost::random::mt19937>(std::time(0));
      this->gen = ::boost::make_shared<boost::random::mt19937>(42);
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
                rdsum+=rd[i]*rd[i];
            }
            rdsum = sqrt(rdsum);
            for(int i=0; i<dim; i++)
            {
                rd[i] /= rdsum;
            }
        }


    bool PcRrt::solve()
    {
      // add the start configuation
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(*this->start));
      // set initial state based on one particle (covariance is zero)
      Particle p(*(this->start));
      ::std::vector<Particle> v;
      v.push_back(p);
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

        ::std::vector<Particle> particles;

        // randomly decide to do a slide or not
        boost::random::uniform_int_distribution<> doSlideDistr(0, 100);
        bool doGuardedMove = doSlideDistr(*this->gen) < 30;
        bool doSlide = doSlideDistr(*this->gen) > 70;
        // sample[...]Particles will return false if the particle set is not useful
        bool sampleResult = false;
        if (doSlide && this->tree[0][n.first].gState->isInCollision())
        {
          sampleResult = this->sampleSlidingParticles(n, chosenSample, this->nrParticles, particles);
        }
        else if (doGuardedMove)
        {
          sampleResult = this->sampleGuardedParticles(n, chosenSample, this->nrParticles, particles);
        }
        else
        {
          sampleResult = this->sampleConnectParticles(n, chosenSample, this->nrParticles, particles);
        }


        if (sampleResult)
        {
          // visualize particles
          this->drawParticles(particles);

          ::boost::shared_ptr<GaussianState> gaussianState = ::boost::make_shared<GaussianState>(particles);
          Gaussian g = gaussianState->configGaussian();
          this->drawEigenvectors(g, 1.0);

          // add a new vertex and edge
          VectorPtr mean = ::boost::make_shared<::rl::math::Vector>(gaussianState->configMean());
          Vertex newVertex = this->addVertex(this->tree[0], mean);
          this->model->setPosition(*mean);
          ::boost::shared_ptr<::rl::math::Transform> tptr = ::boost::make_shared<::rl::math::Transform>(this->model->forwardPosition());
          this->tree[0][newVertex].t = tptr;
          this->tree[0][newVertex].gState = gaussianState;

          this->addEdge(chosenVertex, newVertex, this->tree[0]);

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
            if (this->sampleConnectParticles(nearest, *this->goal, this->nrParticles, goalParticles))
            {
              ::boost::shared_ptr<GaussianState> goalState = ::boost::make_shared<GaussianState>(goalParticles);
              ::rl::math::Real error = goalState->configGaussian().eigenvalues().maxCoeff();
              std::cout << "reached goal with error: " << error << " (max allowed: " << this->goalEpsilon << ")" << std::endl;
              if (error < this->goalEpsilon)
              {
                // visualize goal connect step
                this->drawParticles(goalParticles);
                Gaussian g = goalState->configGaussian();
                this->drawEigenvectors(g, 1.0);
                // add goal connect step to tree
                Vertex connected = this->addVertex(this->tree[0], possibleGoal.q);
                this->tree[0][connected].gState = goalState;
                this->addEdge(possibleGoal.neighbor.first, connected, this->tree[0]);
                this->end[0] = connected;
                return true;
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
    void PcRrt::getPath(VectorList& path)
    {
      Vertex i = this->end[0];

      while (i != this->begin[0])
      {      
        ::boost::shared_ptr<GaussianState> gs = this->tree[0][i].gState;
        Particle p = *(gs->getParticles().begin());
        path.push_front(p.config);
        i = ::boost::source(*::boost::in_edges(i, this->tree[0]).first, this->tree[0]);
      }

      ::boost::shared_ptr<GaussianState> gs = this->tree[0][i].gState;
      Particle p = *(gs->getParticles().begin());
      path.push_front(p.config);

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
        //          noisyFrom[1] += stepY;

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
    bool PcRrt::sampleConnectParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles)
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

        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init);

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

          if(!this->solidScene->getCollisionSurfaceNormal(initPoint,normal))
            return false;

          if(collMap.size()>1)
          {
            std::cout<<"Connect move ends in a state with two collisions - this should not happen!"<<std::endl;
            return false;
          }

          std::string s1=collMap.begin()->first.first;
          std::string s2=collMap.begin()->first.second;

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
          rl::math::Vector3 contactPoint = collMap.begin()->second;
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
#define RANDOM_DIRECTION
#ifdef RANDOM_DIRECTION
      ::rl::math::Vector dir(this->model->getDof());
      sampleDirection(dir);
#endif

      while (pIdx < nrParticles)
      {
        rl::math::Vector init(this->model->getDof());
        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init);

        // sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->configMean();

        ::rl::math::Vector initialError = init - mean;

#ifdef RANDOM_DIRECTION
        ::rl::math::Vector target = init + dir + initialError;
#else
        ::rl::math::Vector target = chosen + initialError;

#endif

        ::rl::math::Vector nextStep = init;

        // move until collision
        int steps = 0;
        bool collision = false;
        bool inInitialCollision = true;

        ::rl::sg::solid::Scene::CollisionMap allColls, initColls;

        ::rl::math::Real step = delta;
        ::rl::math::Real distance = this->model->distance(init, target);
        do
        {
          this->model->interpolateNoisy(init, target,  step / distance, motionNoise, nextStep);

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
          rl::math::Vector3 normal;
          if(!this->solidScene->getCollisionSurfaceNormal(initPoint,normal))
            return false;


          if(collMap.size()>1)
          {
            std::cout<<"Guarded move ends in a state with two collisions - this should not happen!"<<std::endl;
            return false;
          }

          std::string s1 = collMap.begin()->first.first;
          std::string s2 = collMap.begin()->first.second;

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
          rl::math::Vector3 contactPoint = collMap.begin()->second;
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

    bool PcRrt::moveConfigOnSurface(const ::rl::math::Vector& config, const ::rl::math::Vector3& pointOnSurface, const ::rl::math::Vector3& normal, ::rl::math::Vector& out)
    {
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
      while(fabs(dist) > 0.01 && steps++<3)
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

      //Move a bit into the surface until we are colliding
      steps = 0;

      bool collision = this->model->isColliding();
      while(!collision && steps++<10)
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
        collision = this->model->isColliding();
      }

      return collision && fabs(dist < 0.01)  ;
    }

    /**
      Samples a set of particles for a sliding move along a surface.
    */
    bool PcRrt::sampleSlidingParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::std::vector<Particle>& particles)
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


      Contact slidingContact = oldParticle.contacts[slidingIdx];
      std::pair<std::string, std::string> slidingPair;
      slidingPair.first  = slidingContact.shape_robot;
      slidingPair.second = slidingContact.shape_env;

      //Transform EE frame to contact point
      this->model->setPosition(oldParticle.config);
      this->model->updateFrames();
      ::rl::math::Transform ee = this->model->forwardPosition();
      ::rl::math::Transform cp = ee;
      cp.translation() = slidingContact.point - ee.translation();
      ::rl::math::Transform cp_neg = ee;
      cp_neg.translation() = ee.translation() -slidingContact.point;
      this->model->updateTool(cp);
      this->model->updateFrames();
      this->model->updateJacobian();

      ::rl::math::Vector3 slidingNormal = slidingContact.normal_env;
      //this->drawSurfaceNormal(slidingContact.point,slidingContact.normal_env);

      this->model->setPosition(chosen);
      this->model->updateFrames();
      //The pose of a random sample
      ::rl::math::Vector3 chosenForwardPos = this->model->forwardPosition().translation();

      //Project it on the sliding surface
      ::rl::math::Vector3 chosen_proj;
      double distToSurface = projectOnSurface(chosenForwardPos, slidingContact.point, slidingNormal, chosen_proj);


      //This is the target direction of our slide
      ::rl::math::Transform goal = *(this->tree[0][nearest.first].t);
      goal.translation() = chosen_proj;
      //goal.rotation() = //leave constant for now!


      rl::sg::solid::Scene::CollisionMap allColls, finalCollisions;

      //The sliding velocity -- same for all particles
      ::rl::math::Vector6 tdot;

      int pIdx = 0;
      while (pIdx < nrParticles)
      {
        // sample a starting point
        ::rl::math::Vector init(this->model->getDof());
        this->tree[0][nearest.first].gState->sampleConfigurationFromParticle(init);

        this->model->setPosition(init);
        this->model->updateFrames();
        if(!this->model->isColliding())
        {
          std::cout<<"no contact after projection";
        }

        if(pIdx == 0)
        {
          ::rl::math::Transform fp = this->model->forwardPosition();
          ::rl::math::transform::toDelta(fp, goal, tdot);
          tdot.normalize();
        }

        //Sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);


        int steps = 0;
        ::rl::math::Vector nextStep = init;

        while (true)
        {
          this->model->setPosition(nextStep);
          this->model->updateFrames();
          this->model->updateJacobian();
          this->model->updateJacobianInverse();


          ::rl::math::Vector qdot(this->model->getDof());
          qdot.setZero();
          this->model->inverseVelocity(tdot, qdot);

          qdot.normalize();
          qdot *= this->delta;

          ::rl::math::Vector newStep(this->model->getDof());
          this->model->interpolateNoisy(nextStep,nextStep+qdot,1,motionNoise,newStep);

          if(!moveConfigOnSurface(newStep, slidingContact.point, slidingContact.normal_env, nextStep))
          {
            std::cout<<"Could not project cfg on surface!!"<<std::endl;
            this->model->updateTool(cp_neg);
            return false;
          }

          //this->viewer->drawConfiguration(nextStep);

          steps++;

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


        // some magic number
        if (steps < 3)
        {
          // we did not move very far, this is considered failure
          this->model->updateTool(cp_neg);
          return false;
        }

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
            this->model->updateTool(cp_neg);
            return false;
          }
        }

        // valid particle, store it
        if(allColls.size() == 0)
        {
          Particle p(nextStep);
          particles.push_back(p);
        }
        else
        {
          std::vector<Contact> contacts;
          for(::rl::sg::solid::Scene::CollisionMap::iterator it = allColls.begin(); it != allColls.end(); it++)
          {
            std::string c1 = it->first.first;
            std::string c2 = it->first.second;

            ::rl::math::Vector3 contactNormal;
            if(c1==slidingPair.first && c2==slidingPair.second)
            {
              contactNormal = slidingNormal;
            }
            else
            {
              this->solidScene->getCollisionSurfaceNormal(it->second+slidingNormal*this->delta,contactNormal);
            }

            contacts.push_back(Contact(it->second,contactNormal,c1,c2));
          }
          Particle p(nextStep, contacts);
          particles.push_back(p);
        }
        pIdx++;
      }

      this->model->updateTool(cp_neg);
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
        this->model->setPosition(particles[i].config);
        this->model->updateFrames();
        const ::rl::math::Transform& t = this->model->forwardPosition();
        ::rl::math::Vector3 opPos;
        opPos[0] = t.translation().x();
        opPos[1] = t.translation().y();
        opPos[2] = t.translation().z();

        this->viewer->drawSphere(opPos, 0.02);
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
