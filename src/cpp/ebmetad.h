//
// Created by Jennifer Hays on 09/14/2018
//

#ifndef EBMETAD_ENSEMBLEPOTENTIAL_H
#define EBMETAD_ENSEMBLEPOTENTIAL_H

/*! \file
 * \brief Provide EBMetaD MD potential for GROMACS plugin.
 */

#include <vector>
#include <array>
#include <mutex>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/real.h"

// We do not require C++14, so we have a back-ported C++14 feature for C++11 code.
#include "make_unique.h"
#include "sessionresources.h"

// Extra EBMetaD includes
#include <map>

namespace plugin {

    // use the row index for *current* distance, the col index for *historical* distances
    using ForceTable = std::vector<std::vector<double>>;

    struct ebmetad_input {
        /// distance histogram parameters
        double binWidth{0.};

        // historical data
        std::vector<unsigned long> distanceCounts{};
        std::string historicalDataFilename;
        double samplePeriod{0.};

        // This is a table of force magnitudes for faster calculation of the typical MetaDynamics forces.
        // How these force factors should be calculated:
        // 1. Calculate grad(V(collective variable, t))
        // 2. Choose bin size for collective variable.
        // 3. Bin the results of grad(V) according to bin size.
        ForceTable forces{};

        double k{100};
        double maxDist{0.};
        double minDist{0.};

    };

// Generate input structure
    std::unique_ptr<ebmetad_input>
    makeEBMetaDParams(const ForceTable &forces,
                      std::vector<unsigned long> &counts,
                      double binWidth,
                      double samplePeriod,
                      std::string historicalDataFilename,
                      double minDist,
                      double maxDist,
                      double k);

/*!
 * \brief a residue-pair bias calculator for use in EBMetaD ensemble simulations.
 *
 * Applies a force between two sites according to the difference between an experimentally observed
 * site pair distance distribution and the distance distribution observed earlier in the simulation
 * trajectory. This is calculated according to doi: 10.1016/j.bpj.2015.05.024.
 *
 */
    class EBMetaD {
    public:
        using input_param_type = ebmetad_input;

        /* No default constructor. Parameters must be provided. */

        explicit EBMetaD(const input_param_type &params);

        /*!
         * \brief Deprecated constructor taking a parameter list.
         */
        EBMetaD(ForceTable forces,
                std::vector<unsigned long> counts,
                double binWidth,
                double samplePeriod,
                std::string historicalDataFilename,
                double minDist,
                double maxDist,
                double k);

        /*!
         * \brief Evaluates the pair restraint potential.
         *
         * In parallel simulations, the gmxapi framework does not make guarantees about where or
         * how many times this function is called. It should be simple and stateless; it should not
         * update class member data (see ``ensemblepotential.cpp``. For a more controlled API hook
         * and to manage state in the object, use ``callback()``.
         *
         * \param v position of the site for which force is being calculated.
         * \param v0 reference site (other member of the pair).
         * \param t current simulation time (ps).
         * \return container for force and potential energy data.
         */

        gmx::PotentialPointData calculate(gmx::Vector v,
                                          gmx::Vector v0,
                                          gmx_unused double t);

        /*!
         * \brief An update function to be called on the simulation master rank/thread periodically by the Restraint framework.
         *
         * Defining this function in a plugin potential is optional. If the function is defined,
         * the restraint framework calls this function (on the first rank only in a parallel simulation) before calling calculate().
         *
         * The callback may use resources provided by the Session in the callback to perform updates
         * to the local or global state of an ensemble of simulations. Future gmxapi releases will
         * include additional optimizations, allowing call-back frequency to be expressed, and more
         * general Session resources, as well as more flexible call signatures.
         */
        void callback(gmx::Vector v,
                      gmx::Vector v0,
                      double t,
                      const EnsembleResources &resources);

    private:
        /// Width of bins (distance) in histogram
        double binWidth_;
        /// A vector of counts of the number of time that particular distance (bin) has been observed in simulation.
        std::vector<unsigned long> distanceCounts_;
        /// Path to file containing the historical data (distanceCounts_).
        std::string historicalDataFilename_;
        /// Table of force magnitudes
        ForceTable forces_;
        /// Period (in ps) for updating the distanceCounts_ vector. This determines how often to "drop a Gaussian" as in
        /// the standard MetaDynamics procedure.
        double samplePeriod_;

        /// Define the "boundary conditions" for regions outside of the DEER distribution.
        /// Will be used to apply linear biasing potential to keep simulation probability mass inside the DEER region.
        double minDist_;
        double maxDist_;
        double k_;
    };

/*!
 * \brief Use EBMetaD to implement a RestraintPotential
 *
 * This is boiler plate that will be templated and moved.
 */
    class EBMetaDRestraint : public ::gmx::IRestraintPotential, private EBMetaD {
    public:
        using EBMetaD::input_param_type;

        EBMetaDRestraint(std::vector<unsigned long> sites,
                          const input_param_type &params,
                          std::shared_ptr<EnsembleResources> resources
        ) :
                EBMetaD(params),
                sites_{std::move(sites)},
                resources_{std::move(resources)} {}

        /*!
         * \brief Implement required interface of gmx::IRestraintPotential
         *
         * \return list of configured site indices.
         *
         * \todo remove to template header
         * \todo abstraction of site references
         */
        std::vector<unsigned long int> sites() const override {
            return sites_;
        }

        /*!
         * \brief Implement the interface gmx::IRestraintPotential
         *
         * Dispatch to calculate() method.
         *
         * \param r1 coordinate of first site
         * \param r2 reference coordinate (second site)
         * \param t simulation time
         * \return calculated force and energy
         *
         * \todo remove to template header.
         */
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override {
            return calculate(r1,
                             r2,
                             t);
        };

        /*!
         * \brief An update function to be called on the simulation master rank/thread periodically by the Restraint framework.
         *
         * Implements optional override of gmx::IRestraintPotential::update
         *
         * This boilerplate will disappear into the Restraint template in an upcoming gmxapi release.
         */
        void update(gmx::Vector v,
                    gmx::Vector v0,
                    double t) override {
            // Todo: use a callback period to mostly bypass this and avoid excessive mutex locking.
            callback(v,
                     v0,
                     t,
                     *resources_);
        };

        /*!
         * \brief Implement the binding protocol that allows access to Session resources.
         *
         * The client receives a non-owning pointer to the session and cannot extent the life of the session. In
         * the future we can use a more formal handle mechanism.
         *
         * \param session pointer to the current session
         */
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


// Important: Just declare the template instantiation here for client code.
// We will explicitly instantiate a definition in the .cpp file where the input_param_type is defined.
    extern template
    class RestraintModule<EBMetaDRestraint>;

} // end namespace plugin

#endif //EBMETAD_ENSEMBLEPOTENTIAL_H
