//
// Created by Jennifer Hays on 6/12/18.
//

#ifndef GROMACS_LINEARPOTENTIAL_H
#define GROMACS_LINEARPOTENTIAL_H

#include <array>
#include <iostream>
#include <mutex>
#include <vector>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/session.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

#include "make_unique.h"
#include "sessionresources.h"

namespace plugin {

    struct linear_input_param_type {
        double alpha{0};
        double target{0};
        double samplePeriod{0};
        std::string logging_filename{""};
    };

    std::unique_ptr<linear_input_param_type>
    makeLinearParams(double alpha, double target, double samplePeriod, std::string logging_filename);

    class Linear {
    public:
        using input_param_type = linear_input_param_type;

        explicit Linear(const input_param_type &params);

        Linear(double alpha, double target, double samplePeriod,
               std::string filename);

        gmx::PotentialPointData calculate(gmx::Vector v, gmx::Vector v0,
                                          gmx_unused double t);

        void writeparameters(double t, const double R);

        void callback(gmx::Vector v, gmx::Vector v0, double t,
                      const EnsembleResources &resources);

        virtual ~Linear() {}

    private:
        bool initialized_{FALSE};
        double alpha_;

        /// target distance
        double target_;

        // Sample interval
        double samplePeriod_;
        double startTime_{0};
        double nextSampleTime_{0};
        unsigned int currentSample_{0};

        std::string logging_filename_;
        std::unique_ptr<RAIIFile> logging_file_{nullptr};
    };

    class LinearRestraint : public ::gmx::IRestraintPotential,
                            private Linear {
    public:
        using Linear::input_param_type;

        LinearRestraint(const std::vector<unsigned long> &sites,
                        const input_param_type &params,
                        std::shared_ptr<EnsembleResources> resources)
                : Linear(params), sites_{sites}, resources_{std::move(resources)} {}

        std::vector<unsigned long int> sites() const override { return sites_; }

        gmx::PotentialPointData evaluate(gmx::Vector r1, gmx::Vector r2,
                                         double t) override {
            return calculate(r1, r2, t);
        };

        void update(gmx::Vector v, gmx::Vector v0, double t) override {
            callback(v, v0, t, *resources_);
        };

        void bindSession(gmxapi::SessionResources *session) override {
            resources_->setSession(session);
        }

        void setResources(std::unique_ptr<EnsembleResources> &&resources) {
            resources_ = std::move(resources);
        }


    private:
        std::vector<unsigned long int> sites_;
        std::shared_ptr<EnsembleResources> resources_;
    };

    extern template class RestraintModule<LinearRestraint>;
}

#endif // GROMACS_LINEARPOTENTIAL_H

////
//// Created by Jennifer Hays on 5/3/2018
////
//
//#ifndef GROMACS_LINEARPOTENTIAL_H
//#define GROMACS_LINEARPOTENTIAL_H
//
//#include <iostream>
//
//#include "gmxapi/gromacsfwd.h"
//#include "gmxapi/md/mdmodule.h"
//
//#include "gromacs/restraint/restraintpotential.h"
//#include "gromacs/utility/real.h"
//
//namespace plugin {
//
//class Linear {
//public:
//  Linear(real equilibrium, real couplingconstant, std::string logging_filename)
//      : R0{equilibrium}, k{couplingconstant} , logging_filename{logging_filename};
//
//  Linear() : Linear{0.0, 0.0, ""} {};
//
//  // Allow easier automatic generation of bindings.
//  struct input_param_type {
//    float whateverIwant;
//  };
//
//  struct output_type {};
//
//  /*!
//   * \brief Calculate Linear force on particle at position v in reference to
//   * position v0. U = k |v-v0| \param v position at which to evaluate force
//   * \param v0 position of Linear bond reference
//   * \return F = -k (v - v0)/|v - v0|;
//   *
//   * R0 == 1.0 is the equilibrium distance in the Linear potential.
//   * k == 1.0 is the spring constant.
//   *
//   * In the case of a pair of Linearally bonded particles, the force on particle
//   * i is evaluated with particle j as the reference point with \code auto force
//   * = calculateForce(r_i, r_j); \endcode
//   *
//   * The force on particle j is the opposite as the force vector for particle i.
//   * E.g. \code assert(-1 * force, calculateForce(r_j, r_i)); \endcode
//   */
//  gmx::PotentialPointData calculate(gmx::Vector v, gmx::Vector v0,
//                                    gmx_unused double t);
//
//  // Cache of historical distance data. Not thread safe
//  //        std::vector<float> history{};
//
//  // The class will either be inherited as a mix-in or inherit a CRTP base
//  // class. Either way, it probably needs proper virtual destructor management.
//  virtual ~Linear() {
//    //            for (auto&& distance: history)
//    //            {
//    //                std::cout << distance << "\n";
//    //            }
//    //            std::cout << std::endl;
//  }
//
//private:
//  // set equilibrium separation distance
//  // TODO: be clearer about units
//  real R0;
//  // set spring constant
//  // TODO: be clearer about units
//  real k;
//};
//
//// implement IRestraintPotential in terms of Linear
//// To be templated and moved.
//class LinearRestraint : public ::gmx::IRestraintPotential, private Linear {
//public:
//  LinearRestraint(const std::vector<unsigned long> &sites, real R0, real k, std::string logging_filename)
//      : Linear{R0, k, logging_filename}, sites_{sites} {};
//
//  std::vector<unsigned long int> sites() const override;
//
//  // \todo provide this facility automatically
//  gmx::PotentialPointData evaluate(gmx::Vector r1, gmx::Vector r2,
//                                   double t) override;
//
//private:
//  std::vector<unsigned long> sites_;
//};
//
//class LinearModule : public gmxapi::MDModule {
//public:
//  using param_t = Linear::input_param_type;
//
//  LinearModule(std::vector<unsigned long int> sites, real R0, real k, std::string logging_filename) {
//    sites_ = sites;
//    R0_ = R0;
//    k_ = k;
//    logging_filename_ = logging_filename;
//  }
//
//  const char *name() override { return "LinearModule"; }
//
//  /*!
//   * \brief implement gmxapi::MDModule::getRestraint()
//   *
//   * \return Handle to configured library object.
//   */
//  std::shared_ptr<gmx::IRestraintPotential> getRestraint() override {
//    auto restraint = std::make_shared<LinearRestraint>(sites_, R0_, k_);
//    return restraint;
//  }
//
//  /*!
//   * \brief Set restraint parameters.
//   *
//   * \todo generalize this
//   * \param site1
//   * \param site2
//   * \param k
//   * \param R0
//   */
//  void setParams(std::vector<unsigned long int> sites, real R0, real k, std::string logging_filename) {
//    sites_ = sites;
//    R0_ = R0;
//    k_ = k;
//    logging_filename_ = logging_filename;
//  }
//
//private:
//  std::vector<unsigned long int> sites_;
//  real R0_;
//  real k_;
//  std::string logging_filename_;
//};
//
//} // end namespace plugin
//
//#endif // GROMACS_LINEARPOTENTIAL_H
