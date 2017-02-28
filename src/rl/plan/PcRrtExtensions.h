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

      Gaussian()
      {
      }

      Gaussian(const ::std::vector<::rl::math::Vector>& values)
      {
        this->init(values);
      }

      Gaussian(const ::rl::math::Matrix& values)
      {
        this->init(values);
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

      void init(const ::std::vector<::rl::math::Vector>& values)
      {
        ::rl::math::Matrix p;
        p.resize(values.size(), values[0].size());
        for (int i = 0; i < values.size(); ++i)
        {
          p.row(i) = values[i];
        }
        this->init(p);
      }

      void init(const ::rl::math::Matrix& values)
      {
        this->mean = values.colwise().mean();
        if (values.rows() == 1)
        {
          // just one particle
          this->covariance = ::rl::math::Matrix::Zero(values.cols(), values.cols());
          return;
        }
        // substract mean
        ::rl::math::Matrix centered = values.rowwise() - this->mean.transpose();
        // calculate sample covariance
        this->covariance = centered.transpose() * centered / (values.rows()-1);
      }
    };

    class Contact
    {
    public:
      Contact(){}

      Contact(::rl::math::Vector3& p, ::rl::math::Vector3& n, std::string& s_robot, std::string& s_env) :
        point(p),
        normal_env(n),
        shape_robot(s_robot),
        shape_env(s_env)
      {}

      ::rl::math::Vector3 point;
      std::string shape_robot;
      std::string shape_env;
      ::rl::math::Vector3 normal_env; //Normal vector on environment
    };

    class Particle
    {
    public:
      Particle(::rl::math::Vector& q, std::vector<Contact> &c) :
        config(q),
        contacts(c)
      {
      }

    public:
      Particle(::rl::math::Vector& q) :
        config(q)
      {
      }

      ::rl::math::Vector config;
      std::vector<Contact> contacts;
    };

    class GaussianState
    {
    public:
      GaussianState(){}

      GaussianState(const ::std::vector<Particle>& particles, const ::rl::math::Vector3 normal = ::rl::math::Vector3()) :
      gen(42),
      particles(particles),
      slidingNormal(normal)
      {
        this->init(particles);
      }

      void sampleConfiguration(::rl::math::Vector& q)
      {
        assert(q.size() == this->dims);

        for (int i = 0; i < this->dims; ++i)
        {
          q[i] = this->configSampleDistributions[i](this->gen);
        }
      }


      void sampleConfigurationFromParticle(::rl::math::Vector& q, int id)
      {
        assert(q.size() == this->dims);
        //boost::random::uniform_int_distribution<> particleDistr(0, particles.size()-1);
        q = this->particles[id%particles.size()].config;
      }


      Gaussian configGaussian()
      {
        return this->configDistr;
      }

      ::rl::math::Vector configMean()
      {
        return this->configDistr.mean;
      }

      ::rl::math::Matrix configCovariance()
      {
        return this->configDistr.covariance;
      }

      bool isInCollision()
      {
        return (this->particles[0].contacts.size() != 0);
      }

      bool isSlidingMove()
      {
        return (this->slidingNormal.norm()>0.1);
      }

      ::rl::math::Vector3 getSlidingNormal()
      {
        return this->slidingNormal;
      }

      const ::std::vector<Particle>& getParticles()
      {
        return particles;
      }

    private:
      void init(const ::std::vector<Particle>& particles)
      {
        std::vector<::rl::math::Vector> q;
        for(int i=0; i<particles.size(); i++)
        {
          q.push_back(particles[i].config);
        }
        configDistr.init(q);

        this->dims = this->configDistr.covariance.rows();
        for (int i = 0; i < this->dims; ++i)
        {
          ::rl::math::Real mean = this->configDistr.mean[i];
          ::rl::math::Real std_dev = sqrt(this->configDistr.covariance(i,i));
          configSampleDistributions.push_back(boost::random::normal_distribution<>(mean, std_dev));
        }

      }

      ::std::vector<Particle> particles;
      ::rl::math::Vector3 slidingNormal;
      Gaussian configDistr;
      boost::random::mt19937 gen;
      ::std::vector<boost::random::normal_distribution<> > configSampleDistributions;
      int dims;
    };
  }
}

#endif // _RL_PLAN_PCRRTEXTENSIONS_H_
