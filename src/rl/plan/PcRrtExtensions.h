#ifndef _RL_PLAN_PCRRTEXTENSIONS_H_
#define _RL_PLAN_PCRRTEXTENSIONS_H_

#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>

namespace rl
{
  namespace plan
  {
    class Gaussian {
    public:
      Gaussian(const ::rl::math::Matrix& particles)
      {
        this->init(particles);
      }

      Gaussian(const ::std::vector<::rl::math::Vector>& particles)
      {
        ::rl::math::Matrix p;
        p.resize(particles.size(), particles[0].size());
        for (int i = 0; i < particles.size(); ++i)
        {
          p.row(i) = particles[i];
        }
        this->init(p);
      }

      bool isUnimodal()
      {
        return true;
      }

      ::rl::math::Real likelihood(const ::std::vector<::rl::math::Vector>& data)
      {
        ::rl::math::Real det = this->covariance.determinant();
        assert(det != 0.0);

        ::rl::math::Real lh = 0.0;
        for (const auto& d : data)
        {
          ::rl::math::Real p = exp(-0.5*(d - this->mean).transpose()*this->covariance.inverse()*(d - this->mean)) / sqrt((2*M_PI*this->covariance).determinant());
          lh += log(p);
        }

        return lh;
      }

      ::rl::math::Vector mean;
      ::rl::math::Matrix covariance;

    private:
      void init(const ::rl::math::Matrix& particles)
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
    };

    class GaussianState
    {
    public:
      GaussianState(const ::rl::math::Matrix& particles) :
      gaussian(particles),
      gen(42)
      {
        this->init();
      }

      GaussianState(const ::std::vector<::rl::math::Vector>& particles) :
      gaussian(particles),
      gen(42)
      {
        this->init();
      }

      void sample(::rl::math::Vector& q)
      {
        assert(q.size() == this->dims);

        for (int i = 0; i < this->dims; ++i)
        {
          q[i] = this->distributions[i](this->gen);
        }
      }

      ::rl::math::Vector mean()
      {
        return this->gaussian.mean;
      }

      ::rl::math::Matrix covariance()
      {
        return this->gaussian.covariance;
      }
    private:
      void init()
      {
        this->dims = this->gaussian.covariance.rows();
        for (int i = 0; i < this->dims; ++i)
        {
          ::rl::math::Real mean = this->gaussian.mean[i];
          ::rl::math::Real variance = this->gaussian.covariance(i,i);
          distributions.push_back(boost::random::normal_distribution<>(mean, variance));
        }
      }

      Gaussian gaussian;
      boost::random::mt19937 gen;
      ::std::vector<boost::random::normal_distribution<> > distributions;
      int dims;
    };
  }
}

#endif // _RL_PLAN_PCRRTEXTENSIONS_H_