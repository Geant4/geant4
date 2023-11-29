// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Helper namespace field_utils implementation
//
// Author: Dmitry Sorokin, Google Summer of Code 2017
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------

#include "G4FieldUtils.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

namespace field_utils {

G4double absoluteError(const G4double y[], 
                       const G4double yError[],
                             G4double hstep)
{
   const G4double momentum2 = getValue2(y, Value3D::Momentum);
   const G4double invMomentum2 = 1.0 / momentum2; 
   const G4double positionError2 = getValue2(yError, Value3D::Position);
   const G4double momentumError2 = getValue2(yError, Value3D::Momentum);
   const G4double relativeMomentumError2 = momentumError2 * invMomentum2;

   return std::max(std::sqrt(positionError2),
                   std::sqrt(relativeMomentumError2) * hstep);
}

G4double relativeError2(const G4double y[],
                        const G4double yerr[],
                              G4double h,
                              G4double eps_rel_max)
{
   G4double errmax_sq;

   G4double inv_eps_vel_sq = 1.0 / (eps_rel_max * eps_rel_max);
   G4double errvel_sq = 0.0;    // square of momentum vector difference

   G4double eps_pos = eps_rel_max * h; 
   G4double inv_eps_pos_sq = 1.0 / (eps_pos * eps_pos); 

   // Evaluate accuracy
   //
   G4double errpos_sq = getValue2(yerr, Value3D::Position);
   errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

   // Accuracy for momentum
   //
   G4double magvel_sq = getValue2(y, Value3D::Momentum);
   G4double sumerr_sq = getValue2(yerr, Value3D::Momentum); 
   if (magvel_sq > 0.0)
   {
      errvel_sq = sumerr_sq / magvel_sq; 
   }
   else
   {
      G4Exception("field_utils::relativeError","Field001",
                  JustWarning, "found case of zero momentum");
      errvel_sq = sumerr_sq; 
   }
   errvel_sq *= inv_eps_vel_sq;
   errmax_sq = std::max(errpos_sq, errvel_sq);

   return errmax_sq;
}

G4double relativeError(const G4double y[],
                       const G4double yError[],
                       const G4double h,
                       const G4double errorTolerance)
{
   return std::sqrt(relativeError2(y, yError, h, errorTolerance));
}

void copy(G4double dst[], const G4double src[], std::size_t size)
{
   std::memcpy(dst, src, sizeof(G4double) * size);
}


G4double inverseCurvatureRadius(G4double particleCharge,
                                G4double momentum,
                                G4double BField)
{
   return -c_light * particleCharge * BField / momentum;
}

} // field_utils

