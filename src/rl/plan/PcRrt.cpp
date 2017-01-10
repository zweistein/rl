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
       this->gen = ::boost::make_shared<boost::random::mt19937>(std::time(0));
      //this->gen = ::boost::make_shared<boost::random::mt19937>(42);
    }

    //Needed for Guarded moves
    //    void PcRrt::sampleDirection(::rl::math::Vector& rd)
    //    {
    //        boost::random::uniform_real_distribution<> distr(0, 1);
    //        int dim = rd.rows();
    //        ::rl::math::Vector rd(dim);
    //        double rdsum = 0;
    //        for(int i=0; i<dim; i++)
    //        {
    //            rd[dim] = distr(*this->gen);
    //            rsdum+=rd[dim]*rd[dim];
    //        }
    //        rdsum = sqrt(rdsum);
    //        for(int i=0; i<dim; i++)
    //        {
    //            rd[dim] /= rdsum;
    //        }

    //        return rd;
    //    }


    bool PcRrt::solve()
    {
      // add the start configuation
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(*this->start));
      // set initial state based on one particle (covariance is zero)
      ::std::vector<::rl::math::Vector> v;
      ::std::vector<::rl::math::Vector> c;
      v.push_back(*this->start);
      // uncertainty of 0 at starting configuration
      this->tree[0][this->begin[0]].gState = ::boost::make_shared<GaussianState>(v,c);

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

        ::rl::math::Matrix particles;
        ::rl::math::Matrix collisions;

        // randomly decide to do a slide or not
        boost::random::uniform_int_distribution<> doSlideDistr(0, 100);
        bool doGuardedMove = doSlideDistr(*this->gen) < 10;
        bool doSlide = doSlideDistr(*this->gen) > 70;
        // sample[...]Particles will return false if the particle set is not useful
        bool sampleResult;
        bool isInCollision = false;
        if (doSlide && this->tree[0][n.first].gState->isInCollision())
        {
          sampleResult = this->sampleSlidingParticles(n, chosenSample, this->nrParticles, particles, collisions, isInCollision);
        }
        else if (doGuardedMove)
        {
          sampleResult = this->sampleGuardedParticles(n, chosenSample, this->nrParticles, particles, collisions, isInCollision);
        }
        else
        {
          sampleResult = this->sampleConnectParticles(n, chosenSample, this->nrParticles, particles, collisions, isInCollision);
        }


        if (sampleResult)
        {
          // visualize particles
          this->drawParticles(particles);
          // fit a gaussian to the particles
          Gaussian configGaussian(particles);

          Gaussian collisionPointGaussian;
          if(isInCollision)
          {
            // visualize particles
            //this->drawParticles(collisions);
            // fit a gaussian to the collision points
            collisionPointGaussian = Gaussian(collisions);
            this->drawEigenvectors(collisionPointGaussian);
          }
          else
            this->drawEigenvectors(configGaussian);



          // add a new vertex and edge
          Vertex newVertex = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(configGaussian.mean));
          this->tree[0][newVertex].gState = ::boost::make_shared<GaussianState>(particles, collisions);
          this->tree[0][newVertex].gState->setColliding(isInCollision);
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
            ::rl::math::Matrix goalParticles;
            if (this->sampleGoalParticles(nearest, *this->goal, this->nrParticles, goalParticles))
            {
              Gaussian goalGaussian(goalParticles);
              ::rl::math::Real error = goalGaussian.eigenvalues().maxCoeff();
              std::cout << "reached goal with error: " << error << " (max allowed: " << this->goalEpsilon << ")" << std::endl;
              if (error < this->goalEpsilon)
              {
                // visualize goal connect step
                this->drawParticles(goalParticles);
                this->drawEigenvectors(goalGaussian);
                // add goal connect step to tree
                Vertex connected = this->addVertex(this->tree[0], possibleGoal.q);
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

      this->model->interpolate(*tree[nearest.first].q, chosen, step / distance, *last);
      this->model->setPosition(*last);
      this->model->updateFrames();

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

    /**
            Computes a path from start to goal and incorporates motion error if requested.
          */
    void PcRrt::getPath(VectorList& path)
    {
      Vertex i = this->end[0];

      while (i != this->begin[0])
      {
        path.push_front(*this->tree[0][i].q);
        i = ::boost::source(*::boost::in_edges(i, this->tree[0]).first, this->tree[0]);
      }

      path.push_front(*this->tree[0][i].q);

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

    /**
      Samples a set of particles for a move through free-space to a target configuration (chosen)
    */
    bool PcRrt::sampleConnectParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::rl::math::Matrix& particles, ::rl::math::Matrix& collisions, bool& isInCollision)
    {
      particles.resize(nrParticles, this->model->getDof());
      collisions.resize(nrParticles, 3);

      int rowIdx = 0;

      std::vector<::std::string> shapes1, shapes2;
      std::string shape1 = "";
      std::string shape2 = "";

      while (rowIdx < nrParticles)
      {
        // sample a starting point from the gaussian
        Particle init(this->model->getDof());
        this->tree[0][nearest.first].gState->sample(init);

        Particle nextStep = init;

        // sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->mean();

        ::rl::math::Vector initialError = init - mean;
        ::rl::math::Vector target = chosen + initialError;

        int steps = 0;
        bool reached = false;
        bool collision = false;
        ::rl::math::Real step = this->delta;
        ::rl::math::Real distance = this->model->distance(init, target);
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

          step += delta;
          steps++;
        }
        while (!collision && !reached);

        if (reached && !collision)
        {
          particles.row(rowIdx) = nextStep;
          rowIdx++;
          shape1 = "no collision";
          isInCollision = false;
        }
        else if (steps > 1 && collision)
        {
          std::string s1, s2;
          this->solidScene->lastCollidingShapes(s1, s2);
          if (shape1 == "")
          {
            // first collision, init contact state shapes
            shape1 = s1;
            shape2 = s2;
          }
          else if (s1 != shape1 || s2 != shape2)
          {
            // not every collision is between the same two shapes, -> not one contact state
            return false;
          }

          // valid particle, store it
          particles.row(rowIdx) = nextStep;
          collisions.row(rowIdx) = this->solidScene->lastCollisionPoint();
          rowIdx++;
          isInCollision = true;
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
    bool PcRrt::sampleGuardedParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::rl::math::Matrix& particles, ::rl::math::Matrix& collisions, bool& isInCollision)
    {
      particles.resize(nrParticles, this->model->getDof());
      collisions.resize(nrParticles, 3);

      int rowIdx = 0;

      std::vector<std::string> shapes1, shapes2;
      std::string shape1 = "";
      std::string shape2 = "";

      while (rowIdx < nrParticles)
      {
        // sample a starting point from the gaussian
        Particle init(this->model->getDof());
        this->tree[0][nearest.first].gState->sample(init);

        Particle nextStep = init;

        // sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->mean();

        ::rl::math::Vector initialError = init - mean;
        ::rl::math::Vector target = chosen + initialError;

        // move until collision
        int steps = 0;
        bool collision = false;
        ::rl::math::Real step = delta;
        ::rl::math::Real distance = this->model->distance(init, target);
        do
        {
          this->model->interpolateNoisy(init, target,  step / distance, motionNoise, nextStep);

          this->model->setPosition(nextStep);
          this->model->updateFrames();

          collision = this->model->isColliding();

          step += delta;
          steps++;
        }
        while (!collision);

        if (steps > 1)
        {
          std::string s1, s2;
          this->solidScene->lastCollidingShapes(s1, s2);
          if (shape1 == "")
          {
            // first collision, init contact state shapes
            shape1 = s1;
            shape2 = s2;
          }
          else if (s1 != shape1 || s2 != shape2)
          {
            // not every collision is between the same two shapes, -> not one contact state
            return false;
          }

          // valid particle, store it
          particles.row(rowIdx) = nextStep;
          collisions.row(rowIdx) = this->solidScene->lastCollisionPoint();
          rowIdx++;
        }
        else
        {
          // we moved just one step until the next collision, this is not useful
          return false;
        }
      }

      isInCollision = true;
      return true;
    }

    bool PcRrt::projectOnSurface(const ::rl::math::Vector& point, const ::rl::math::Vector& pointOnSurface, const ::rl::math::Vector& normal, ::rl::math::Vector& out)
    {
      double dist = (point-pointOnSurface).dot(normal);
      out = point-dist*normal;
      return dist;
    }


    /**
      Samples a set of particles for a sliding move along a surface.
    */
    bool PcRrt::sampleSlidingParticles(const Neighbor& nearest, const ::rl::math::Vector& chosen, int nrParticles, ::rl::math::Matrix& particles, ::rl::math::Matrix& collisions, bool& isInCollision)
    {
      if (!this->tree[0][nearest.first].gState->isInCollision())
      {
        // we must start from within a collision to slide
        return false;
      }

      particles.resize(nrParticles, this->model->getDof());
      collisions.resize(nrParticles, 3);

      // project chosen on the plane
      ::rl::math::Vector normal(this->model->getDof());
      if (!this->getNormal(nearest.first, normal))
      {
        // uncertainty distribution insufficient for retrieving normal
        std::cout << "uncertainty distribution insufficient for retrieving normal" << std::endl;
        return false;
      }
      ::rl::math::Vector pVec = *(this->tree[0][nearest.first].q);
      ::rl::math::Vector goal;
      double dist = projectOnSurface(chosen, pVec, normal, goal);

      if (dist < 0)
      {
        dist *= -1;
        normal *= -1;
      }

      // retrieve direction of normal
      ::rl::math::Vector testPoint(this->model->getDof());
      testPoint = pVec + normal*this->delta;
      this->model->setPosition(testPoint);
      this->model->updateFrames();
      if (this->model->isColliding()) 
      {
        // inverse normal
        normal *= -1;
      }
      else
      {
        testPoint = pVec - normal*this->delta;
        this->model->setPosition(testPoint);
        this->model->updateFrames();
        if (!this->model->isColliding())
        {
          // no collision in any direction, this is bad
          std::cout << "error getting normal direction" << std::endl;
          return false;
        }
      }

      //this->drawSurfaceNormal(pVec, normal);

      int rowIdx = 0;
      std::string shape1, shape2;
      std::string finalCollision = "";

      while (rowIdx < nrParticles)
      {
        // sample a starting point
        Particle init(this->model->getDof());

        // must be in collision
        do
        {
          this->tree[0][nearest.first].gState->sample(init);
          projectOnSurface(init, pVec, normal, init);
          this->model->setPosition(init);
          this->model->updateFrames();
        }
        while(!this->model->isColliding());


        //Sample noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        this->model->sampleMotionError(motionNoise);
        projectOnSurface(motionNoise, pVec, normal, motionNoise);

        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->mean();

        std::string slidingSurface = "";
        std::map<::std::string, bool> allColls;

        // aquire sliding surface
        this->getAllCollidingShapes(allColls);
        if (allColls.size() != 1)
        {
          // we should not collide with more than one obstacle here, so this is a weird state
          std::cout << "weird slidingSurface init, abort" << std::endl;
          return false;
        }
        slidingSurface = allColls.begin()->first;

        ::rl::math::Vector initialError = init - mean;
        ::rl::math::Vector target = goal + initialError;

        int steps = 0;
        double step = delta;
        double distance = this->model->distance(init, target);

        ::rl::math::Vector nextStep(this->model->getDof());

        while (true)
        {
          this->model->interpolateNoisy(init, target,  step / distance, motionNoise, nextStep);
          this->model->setPosition(nextStep);
          this->model->updateFrames();

          step += delta;
          steps++;

          if (this->model->isColliding())
          {
            this->getAllCollidingShapes(allColls);

            if (allColls[slidingSurface])
            {
              // we moved into sliding surface, so reflect out
              do
              {
                // move into direction of surface normal to escape from sliding surface
                nextStep += normal * this->delta;
                this->model->setPosition(nextStep);
                this->model->updateFrames();
                this->getAllCollidingShapes(allColls);
              }
              while (allColls[slidingSurface]);

              // check if we are in free-space or if we have another collision
              if (this->model->isColliding())
              {
                // there is still a collision, so the sliding ends here
                break;
              }
            }
            else
            {
              // we are not stuck in the sliding surface, but we are colliding,
              // this move ends here
              break;
            }
          }
          else
          {
            // lost contact
            break;
          }
        }

        // some magic number
        if (steps < 10)
        {
          // we did not move very far, this is considered failure
          return false;
        }

        this->getAllCollidingShapes(allColls);
        if (finalCollision == "")
        {
          // first particle, init final collision
          if (allColls.size() == 0)
          {
            finalCollision = "lost_contact";
          }
          else
          {
            finalCollision = allColls.begin()->first;
          }
        }
        else
        {
          if (allColls.size() == 0)
          {
            if (finalCollision != "lost_contact")
            {
              // some particles lost contact, some didn't...
              return false;
            }
          }
          else
          {
            if (finalCollision != allColls.begin()->first)
            {
              // different collision than the other particles -> different contact state
              return false;
            }
          }
        }

        // valid particle, store it
        particles.row(rowIdx) = nextStep;
        if (finalCollision != "lost_contact")
          collisions.row(rowIdx) = this->solidScene->lastCollisionPoint();
        rowIdx++;
      }

      isInCollision = (finalCollision != "lost_contact");
      return true;
    }

    bool PcRrt::sampleGoalParticles(const Neighbor& nearest, ::rl::math::Vector& goal, int nrParticles, ::rl::math::Matrix& particles)
    {
      double posError = 0.01;
      boost::random::normal_distribution<> distr(0, posError);
      particles.resize(nrParticles, this->model->getDof());

      int rowIdx = 0;

      while (rowIdx < nrParticles)
      {
        Particle nextStep(this->model->getDof());

        // sample initial position
        this->tree[0][nearest.first].gState->sample(nextStep);

        // sample independent noise
        ::rl::math::Vector motionNoise(this->model->getDof());
        ::rl::math::Vector mean = this->tree[0][nearest.first].gState->mean();
        ::rl::math::Vector delta = goal - mean;

        for (int i = 0; i < this->model->getDof(); ++i)
        {
          motionNoise[i] = distr(*this->gen) * delta[i];
        }


        ::rl::math::Vector error = nextStep - mean;
        ::rl::math::Vector target = goal + error + motionNoise;

        int steps = 0;
        bool reached = false;
        bool collision = false;
        do
        {
          ::rl::math::Real distance = this->model->distance(nextStep, target);
          ::rl::math::Real step = distance;

          if (step <= this->delta)
          {
            reached = true;
          }
          else
          {
            step = this->delta;
          }
          this->model->interpolate(nextStep, target,  step / distance, nextStep);

          this->model->setPosition(nextStep);
          this->model->updateFrames();

          collision = this->model->isColliding();
          steps++;
        }
        while (!collision && !reached);

        if (reached)
        {
          particles.row(rowIdx) = nextStep;
          rowIdx++;
        }
        else
        {
          return false;
        }
      }

      return true;
    }

    bool PcRrt::getNormal(const Vertex& vertex, ::rl::math::Vector& normal)
    {
      ::rl::math::Vector evals = this->tree[0][vertex].gState->gaussian().eigenvalues();
      ::rl::math::Matrix evecs = this->tree[0][vertex].gState->gaussian().eigenvectors();

      int minIdx, maxIdx;
      evals.minCoeff(&minIdx);
      evals.maxCoeff(&maxIdx);

      normal = evecs.col(minIdx).normalized();
      // magic number assures that there is sufficient uncertainty in one direction
      return evals[maxIdx] > 0.001;
    }

    void PcRrt::getAllCollidingShapes(::std::map<::std::string, bool>& collidingShapes)
    {
      collidingShapes.clear();
      ::rl::sg::Body *sceneBody = this->solidScene->getModel(1)->getBody(0);
      ::rl::sg::Shape *robotShape = this->solidScene->getModel(0)->getBody(2)->getShape(0);

      ::std::string shape1, shape2;
      for (size_t i = 0; i < sceneBody->getNumShapes(); ++i)
      {

        if (this->solidScene->areColliding(sceneBody->getShape(i), robotShape))
        {
          this->solidScene->lastCollidingShapes(shape1, shape2);
          collidingShapes[shape1] = true;
        }

      }
    }

    void PcRrt::drawParticles(::rl::math::Matrix& particles)
    {
      for (int rowIdx = 0; rowIdx < particles.rows(); ++rowIdx)
      {
        this->model->setPosition(particles.row(rowIdx));
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
      ::rl::math::Vector3 start, end1, end2;
      ::rl::math::Matrix vectors = gaussian.eigenvectors();
      ::rl::math::Vector values = gaussian.eigenvalues();
      ::rl::math::Vector eigen1 = vectors.col(0) * sqrt(values(0)) * scale;
      ::rl::math::Vector eigen2 = vectors.col(1) * sqrt(values(1)) * scale;

      if(gaussian.mean.rows()==3)
      {
        start = gaussian.mean;
        end1 = gaussian.mean + eigen1;
        end2 = gaussian.mean + eigen2;
      }
      else
      {
        this->model->setPosition(gaussian.mean);
        this->model->updateFrames();
        const ::rl::math::Transform& t1 = this->model->forwardPosition();
        start(0) = t1.translation().x();
        start(1) = t1.translation().y();
        start(2) = t1.translation().z();


        this->model->setPosition(gaussian.mean + eigen1);
        this->model->updateFrames();
        const ::rl::math::Transform& t2 = this->model->forwardPosition();
        end1(0) = t2.translation().x();
        end1(1) = t2.translation().y();
        end1(2) = t2.translation().z();


        this->model->setPosition(gaussian.mean + eigen2);
        this->model->updateFrames();
        const ::rl::math::Transform& t3 = this->model->forwardPosition();
        end2(0) = t3.translation().x();
        end2(1) = t3.translation().y();
        end2(2) = t3.translation().z();
      }




      this->viewer->drawLine(start, end1);
      this->viewer->drawLine(start, end2);
    }

    void PcRrt::drawSurfaceNormal(::rl::math::Vector& startPoint, ::rl::math::Vector& normal, ::rl::math::Real scale)
    {
      ::rl::math::Vector3 start, end;

      this->model->setPosition(startPoint);
      this->model->updateFrames();
      const ::rl::math::Transform& t1 = this->model->forwardPosition();
      start(0) = t1.translation().x();
      start(1) = t1.translation().y();
      start(2) = 0.6;

      this->model->setPosition((startPoint + normal) * scale);
      this->model->updateFrames();
      const ::rl::math::Transform& t2 = this->model->forwardPosition();
      end(0) = t2.translation().x();
      end(1) = t2.translation().y();
      end(2) = 0.6;

      this->viewer->drawLine(start, end);
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
