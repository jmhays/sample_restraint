//
// Created by Jennifer Hays on 5/3/2018
//

#include "linearpotential.h"
#include <cmath>

#include <array>

namespace plugin {

gmx::PotentialPointData Linear::calculate(gmx::Vector v, gmx::Vector v0,
                                          gmx_unused double t) {
  auto rdiff = v - v0;
  const auto Rsquared = dot(rdiff, rdiff);
  const auto R = sqrt(Rsquared);
  const auto Rdiff = R - R0;
  const auto Rdiffsqared = Rdiff * Rdiff;
  const auto diffR = sqrt(Rdiffsqared);
  // TODO: find appropriate math header and namespace

  gmx::PotentialPointData output;

  output.energy = k * diffR;
  // Direction of force is ill-defined when v == v0
  if (R != 0 && R != R0) {
    if (R > R0)
      output.force = -(k / R0 / double(R)) * rdiff;
    else
      output.force = (k / R0 / double(R)) * rdiff;
  }
  return output;
}

gmx::PotentialPointData LinearRestraint::evaluate(gmx::Vector r1,
                                                  gmx::Vector r2, double t) {
  return calculate(r1, r2, t);
}

std::vector<unsigned long int> LinearRestraint::sites() const { return sites_; }

} // end namespace plugin
