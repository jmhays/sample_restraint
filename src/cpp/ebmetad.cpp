//
// Created by Jennifer Hays on 09/14/2018
//

/*! \file
 * \brief Code to implement the potential declared in ebmetad.h
 *
 * This file currently contains boilerplate that will not be necessary in future gmxapi releases, as
 * well as additional code used in implementing the EBMetaD workflow.
 *
 */

#include "ebmetad.h"

#include <cmath>

#include <vector>

#include "gmxapi/context.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdsignals.h"

#include "sessionresources.h"
#include <fstream>

namespace plugin {

    unsigned long distToBin(double distance, double binWidth){
        // Convert a distance (calculated in the callback function) to a bin number.
        return (unsigned long)(distance/binWidth);
    }


    void writeHistoricalData(std::string filename, std::vector<unsigned long> distanceCounts){
        std::ofstream outfile(filename);
        for (auto &count: distanceCounts){
            outfile << count << "\n";
        }
        outfile.close();
    }

    EBMetaD::EBMetaD(ForceTable forces,
                     std::vector<unsigned long> distanceCounts,
                     double binWidth,
                     double samplePeriod,
                     std::string historicalDataFilename,
                     double minDist,
                     double maxDist,
                     double k) :
            forces_{std::move(forces)},
            distanceCounts_{std::move(distanceCounts)},
            binWidth_{binWidth},
            samplePeriod_{samplePeriod},
            historicalDataFilename_{std::move(historicalDataFilename)},
            minDist_{minDist},
            maxDist_{maxDist},
            k_{k} {};

    EBMetaD::EBMetaD(const plugin::EBMetaD::input_param_type &params) :
            EBMetaD(params.forces,
                    params.distanceCounts,
                    params.binWidth,
                    params.samplePeriod,
                    params.historicalDataFilename,
                    params.minDist,
                    params.maxDist,
                    params.k) {};

//
//
// HERE is the (optional) function that updates the state of the restraint periodically.
// It is called before calculate() once per timestep per simulation (on the master rank of
// a parallelized simulation).
//
//
    void EBMetaD::callback(gmx::Vector v,
                                    gmx::Vector v0,
                                    double t,
                                    const EnsembleResources &resources) {
        auto rdiff = v - v0;
        const auto Rsquared = dot(rdiff,
                                  rdiff);
        const auto R = sqrt(Rsquared);

        // Store historical data every samplePeriod steps
	// Time for an update if modulus is zero
	auto time_for_update = ! (bool) std::fmod(t, samplePeriod_);
        if (time_for_update) {
            fprintf(stderr, "UPDATING DISTANCE COUNTS: t=%f d=%f\n", t, R);
            distanceCounts_.at(distToBin(R, binWidth_))++;

            // Just for reporting
            double f{0};
            auto forceVector = forces_[distToBin(R, binWidth_)];
            for (unsigned long i = 0; i < distanceCounts_.size(); i++){
                    f += forceVector[i]*distanceCounts_[i];
                }
            fprintf(stderr, "NEW FORCES WILL BE: %f\n", f);
        };

        // Write the new counts to a file so that simulations may be resumed.
        writeHistoricalData(historicalDataFilename_, distanceCounts_);
    }


//
//
// HERE is the function that does the calculation of the restraint force.
//
//
    gmx::PotentialPointData EBMetaD::calculate(gmx::Vector v,
                                                        gmx::Vector v0,
                                                        double t) {
        // This is not the vector from v to v0. It is the position of a site
        // at v, relative to the origin v0. This is a potentially confusing convention...
        auto rdiff = v - v0;
        const auto Rsquared = dot(rdiff,
                                  rdiff);
        const auto R = sqrt(Rsquared);


        // Compute output
        gmx::PotentialPointData output;
        // Energy not needed right now.
        //    output.energy = 0;

        if (R != 0) // Direction of force is ill-defined when v == v0
        {

            double f{0};

            if (R > maxDist_) {
                // apply a force to reduce R
                f = (k_/norm(rdiff)) * (maxDist_ - R);
            } else if (R < minDist_) {
                // apply a force to increase R
                f = (k_/norm(rdiff)) * (minDist_ - R);
            } else {
                /// Here we do EBMetaD evaluation
                /// The force will be a sum of the force factors associated with the current distance: these factors
                /// depend on 1) the historical data, i.e., the distances already visited in the simulation
                /// (distanceCounts_), and 2) the DEER probability distributions.
                auto forceVector = forces_[distToBin(R, binWidth_)];
                for (unsigned long i = 0; i < distanceCounts_.size(); i++){
                    f += forceVector[i]*distanceCounts_[i];
                }
            }

            output.force = f * rdiff;
        }
        return output;
    }

    std::unique_ptr<ebmetad_input>
     makeEBMetaDParams(const ForceTable &forces,
                      std::vector<unsigned long> &distanceCounts,
                      double binWidth,
                      double samplePeriod,
                      std::string &historicalDataFilename,
                      double minDist,
                      double maxDist,
                      double k) {
        using gmx::compat::make_unique;
        auto params = make_unique<ebmetad_input>();

        params->forces = forces;
        params->distanceCounts = distanceCounts;
        params->binWidth = binWidth;
        params->samplePeriod = samplePeriod;
        params->historicalDataFilename = historicalDataFilename;
        params->minDist = minDist;
        params->maxDist = maxDist;
        params->k = k;

        return params;
    };

// Important: Explicitly instantiate a definition for the templated class declared in ebmetad.h.
// Failing to do this will cause a linker error.
    template
    class ::plugin::RestraintModule<EBMetaDRestraint>;

} // end namespace plugin
