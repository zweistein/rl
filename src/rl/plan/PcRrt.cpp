#include <iostream>
#include <cmath>
#include <ctime>

#include <boost/make_shared.hpp>
#include <boost/random.hpp>

#include <Eigen/Eigenvalues>

#include "SimpleModel.h"
#include "Viewer.h"

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
      this->begin[0] = this->addVertex(this->tree[0], ::boost::make_shared< ::rl::math::Vector >(*this->start));
      // set initial covariance
      this->tree[0][this->begin[0]].covariance = ::rl::math::Matrix::Zero(this->model->getDof(), this->model->getDof());

      
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

        // PcRrt::ParticleSet particles;
        ::rl::math::Matrix particles;
        // check if the extension was successful
        if (this->sampleParticles(chosenVertex, angle, this->nrParticles, particles))
        {
          this->drawParticles(particles);
          // fit a gaussian to the particles
          Gaussian gaussian(particles);

          // TESTING
          // std::cout << "all = matrix(c(";
          // for (int i = 0; i < particles.rows(); ++i)
          // {
          //   std::cout << particles.row(i)[0] << "," << particles.row(i)[1];
          //   if (i < particles.rows()-1) std::cout << ",";
          // }
          // std::cout << "), ncol=2, byrow=TRUE)" << std::endl;
          // for (int k = 1; k <= 2; ++k)
          // {
          //   ::std::vector<::std::vector<::rl::math::Vector> > clusters(k);
          //   this->kMeans(particles, k, clusters);
          //   ::rl::math::Real lh = 0.0;
          //   ::rl::math::Real sse = 0.0;
          //   for (int i = 0; i < k; ++i)
          //   {
          //     if (clusters[i].size() == 0)
          //     {
          //       std::cout << "empty" << std::endl;
          //       continue;
          //     } 
          //     Gaussian g(clusters[i]);
          //     lh += g.likelihood(clusters[i]);

          //     // SSE
          //     for (const auto& d : clusters[i])
          //     {
          //       sse += pow(this->model->transformedDistance(d, g.mean), 2);
          //     }
          //   }
            // lh = exp(lh);
          //   std::cout << k << " clusters, lh = " << lh << ", sse = " << sse << std::endl;
          // }


          // for (const auto& c : clusters)
          // {
          //   std::cout << "" << std::endl;
          //   for (const auto& p : c)
          //   {
          //     std::cout << p << std::endl;
          //   }
          // }
          // std::cout << std::endl;
          // ::std::vector<::std::vector<::rl::math::Vector> > clusters(2);
          // this->kMeans(particles, 2, clusters);

          // std::cout << "c1 = matrix(c(";
          // for (int i = 0; i < clusters[0].size(); ++i)
          // {
          //   std::cout << clusters[0][i][0] << "," << clusters[0][i][1];
          //   if (i < clusters[0].size()-1) std::cout << ",";
          // }
          // std::cout << "), ncol=2, byrow=TRUE)" << std::endl;
          // TESTING


          // if (gaussian.isUnimodal())
          // if (clusters[0].size() == 0 || clusters[1].size() == 0)
          if (this->isUnimodal(particles))
          {
            Vertex newVertex = this->addVertex(this->tree[0], ::boost::make_shared<::rl::math::Vector>(gaussian.mean));
            this->tree[0][newVertex].covariance = gaussian.covariance;

            this->addEdge(chosenVertex, newVertex, this->tree[0]);

            // try to connect to goal
            Neighbor nearest;
            nearest.first = newVertex;
            nearest.second = this->model->transformedDistance(*this->tree[0][newVertex].q, *this->goal);

            PossibleGoal possibleGoal;
            possibleGoal.neighbor = nearest;
            possibleGoal.q = this->tryConnect(this->tree[0], nearest, *this->goal);

            if (NULL != possibleGoal.q && this->areEqual(*possibleGoal.q, *this->goal)) {
              Vertex connected = this->addVertex(this->tree[0], possibleGoal.q);
              this->addEdge(possibleGoal.neighbor.first, connected, this->tree[0]);
              this->end[0] = connected;
              return true;
            }
          }

          int foo;
          // std::cin >> foo;
        }



        // if (NULL != extended)
        // {
        //   if (this->areEqual(*this->tree[0][extended].q, *this->goal))
        //   {
        //     this->end[0] = extended;
        //     return true;
        //   }
        // }
        
        timer.stop();
      }
      
      return false;
    }

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
            // std::cout << "k=" << k << " " << i+1 << " empty" << std::endl;
            emptyClusterCount++;
          }
          // for every data point in cluster[i]
          // for (const auto& p : clusters[i])
          // {
          //   Gaussian g(clusters[i]);
          //   sse[k-1] += pow(this->model->transformedDistance(p, g.mean), 2);
          // }
        }
      }

      // ::rl::math::Real sseRatio = (sse[1] - sse[0]) / (sse[2] - sse[1]);
      // std::cout << "sse-ratio = " << sseRatio << std::endl;

      // for (int i = 0; i < maxK; ++i)
      // {
      //   // ::rl::math::Real sseRatio = sse[i-1] / sse[i];
      //   // std::cout << i-1 << "/" << i << " sse-ratio = " << sseRatio << std::endl;
      //   std::cout << "k=" << i+1 << ", sse = " << sse[i] << std::endl;
      // }

      if (emptyClusterCount == maxK*(maxK+1)/2 - maxK) {
        std::cout << "good" << std::endl;
        return true;
      } else {
        std::cout << "bad" << std::endl;
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

      // Vertex connected = this->addVertex(tree, last);
      // this->addEdge(nearest.first, connected, tree);
      // return connected;
    }

    // ::boost::shared_ptr<PcRrt::ParticleSet> PcRrt::sampleParticles(const Vertex& start, float angle, int nrParticles)
    bool PcRrt::sampleParticles(const Vertex& start, float angle, int nrParticles, ::rl::math::Matrix& particles)
    {
      // ::boost::shared_ptr<PcRrt::ParticleSet> particles = ::boost::make_shared<PcRrt::ParticleSet>(nrParticles);

      // boost::random::mt19937 gen(std::time(0));
      boost::random::normal_distribution<> distr(angle, 10*M_PI/360);
      
      particles.resize(nrParticles, this->model->getDof());

      int fails = 0;
      int rowIdx = 0;

      while (rowIdx < nrParticles)
      {
        Particle nextStep(*this->tree[0][start].q);
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
          particles.row(rowIdx) = nextStep;
          rowIdx++;
        }
        else
        {
          fails++;
          if (fails > nrParticles / 4)
          {
            // std::cout << "over" << std::endl;
            return false;
          }
        }
      }
      return true;
      // return particles;
    }
    
    ::rl::math::Vector PcRrt::sampleDirection(Vertex& vertex)
    {
      // calculate eigenvectors
      ::Eigen::EigenSolver<::rl::math::Matrix> eig(this->tree[0][vertex].covariance);
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

        // update means
        for (int i = 0; i < k; ++i)
        {
          // if (clusters[i].empty())
          // {
          //   // if cluster is empty, keep mean from last time
          //   continue;
          // }
          means[i] = ::rl::math::Vector(this->model->getDof());
          for (const auto& p : clusters[i])
          {
            means[i] += p;
          }
          means[i] /= clusters[i].size();
        }
      }
    }
  }
}