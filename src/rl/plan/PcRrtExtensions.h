#ifndef _RL_PLAN_PCRRTEXTENSIONS_H_
#define _RL_PLAN_PCRRTEXTENSIONS_H_

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

      ::rl::math::Real mahalanobis(::rl::math::Vector& x)
      {
        return sqrt((x - this->mean).transpose() * this->covariance.inverse() * (x - this->mean));
      }

      ::rl::math::Matrix eigenvectors()
      {
        ::Eigen::EigenSolver<::rl::math::Matrix> eig(this->covariance);
        return eig.eigenvectors().real();
      }

      ::rl::math::Vector eigenvalues()
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
      inCollision(false)
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

      Gaussian gaussian()
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

      bool isInCollision()
      {
        return this->inCollision;
      }

    private:
      ::rl::math::Matrix part;
      Gaussian gaussianDistr;
      boost::random::mt19937 gen;
      ::std::vector<boost::random::normal_distribution<> > distributions;
      int dims;
      bool inCollision;
    };
  }
}

#endif // _RL_PLAN_PCRRTEXTENSIONS_H_
