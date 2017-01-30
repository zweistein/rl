#ifndef _RL_PLAN_PCRRTEXTENSIONS_H_
#define _RL_PLAN_PCRRTEXTENSIONS_H_

#include <cmath>
#include <set>

#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

#include <Eigen/Eigenvalues>

namespace rl
{
  namespace plan
  {
    class Gaussian {
    public:
      Gaussian(const ::rl::math::Matrix& particles)
      {
        this->mean = particles.colwise().mean();
        if (particles.rows() == 1)
        {
          // just one particle
          this->covariance = ::rl::math::Matrix::Zero(particles.cols(), particles.cols());
          return;
        }
        // substract mean
        ::rl::math::Matrix centered = particles.rowwise() - this->mean.transpose();
        // calculate sample covariance
        this->covariance = centered.transpose() * centered / (particles.rows()-1);
      }

      ::rl::math::Real mahalanobis(::rl::math::Vector& x) const
      {
        return sqrt((x - this->mean).transpose() * this->covariance.inverse() * (x - this->mean));
      }

      ::rl::math::Matrix eigenvectors() const
      {
        ::Eigen::EigenSolver<::rl::math::Matrix> eig(this->covariance);
        return eig.eigenvectors().real();
      }

      ::rl::math::Vector eigenvalues() const
      {
        ::Eigen::EigenSolver<::rl::math::Matrix> eig(this->covariance);
        return eig.eigenvalues().real();
      }

      ::rl::math::Vector mean;
      ::rl::math::Matrix covariance;
    };

    class GaussianState
    {
    public:
      GaussianState(const ::rl::math::Matrix& particles) :
      part(particles),
      gaussianDistr(particles),
      gen(42),
      inCollision(false),
      normalSet(false)
      {
        this->dims = this->gaussianDistr.covariance.rows();
        for (int i = 0; i < this->dims; ++i)
        {
          ::rl::math::Real mean = this->gaussianDistr.mean[i];
          ::rl::math::Real std_dev = sqrt(this->gaussianDistr.covariance(i,i));
          distributions.push_back(boost::random::normal_distribution<>(mean, std_dev));
        }
      }

      void sample(::rl::math::Vector& q)
      {
        assert(q.size() == this->dims);

        for (int i = 0; i < this->dims; ++i)
        {
          q[i] = this->distributions[i](this->gen);
        }
      }

      Gaussian gaussian() const
      {
        return this->gaussianDistr;
      }

      ::rl::math::Matrix& particles()
      {
        return this->part;
      }

      ::rl::math::Vector& mean()
      {
        return this->gaussianDistr.mean;
      }

      ::rl::math::Matrix& covariance()
      {
        return this->gaussianDistr.covariance;
      }

      void setColliding(bool isInCollision)
      {
        this->inCollision = isInCollision;
      }

      bool isInCollision() const
      {
        return this->inCollision;
      }

      ::rl::math::Vector3& getNormal()
      {
        return this->normal;
      }

      void setNormal(::rl::math::Vector3 normal)
      {
        this->normalSet = true;
        this->normal = normal;
      }

      bool hasNormal() const
      {
        return this->normalSet;
      }

    private:
      ::rl::math::Matrix part;
      ::rl::math::Vector3 normal;
      Gaussian gaussianDistr;
      boost::random::mt19937 gen;
      ::std::vector<boost::random::normal_distribution<> > distributions;
      int dims;
      bool inCollision;
      bool normalSet;
    };

    class ExtensionState
    {
    public:
      void setGuardedSlide(::rl::math::Vector& dir)
      {
        int angle = this->vectorToAngle(dir);
        executedDirection.insert(angle);
      }

      bool executedGuardedSlide(::rl::math::Vector& dir)
      {
        return executedDirection.find(this->vectorToAngle(dir)) != executedDirection.end();
      }

    private:
      int vectorToAngle(::rl::math::Vector& vec)
      {
        return (std::atan2(vec[1], vec[0]) - std::atan2(0, 1)) / M_PI * 180 / 45;
      }

      std::set<int> executedDirection;
    };
  }
}

#endif // _RL_PLAN_PCRRTEXTENSIONS_H_
