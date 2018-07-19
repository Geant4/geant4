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
//
// $Id: $
//
//
//    Implementation by Dmitry Sorokin - GSoC 2017
//       Work supported by Google as part of Google Summer of Code 2017.
//    Supervision / code review: John Apostolakis


#include "G4FieldUtils.hh"

namespace field_utils {

G4double relativeError(
    const G4double y[],
    const G4double yError[],
    const G4double h,
    const G4double errorTolerance)
{
    // Accuracy for position
    G4double error2 = getValue2(yError, Value3D::Position) / sqr(h);

    // Accuracy for momentum
    const G4double momentum2 = getValue2(y, Value3D::Momentum);
    if (momentum2 > 0) {
       const G4double momentumError2 =
           getValue2(yError,  Value3D::Momentum) / momentum2;
       error2 = std::max(error2, momentumError2);
    } else {
        G4Exception("field_utils::relativeError","Field001",
                    JustWarning, "found case of zero momentum");
    }
#if 0
    // Accuracy for spin
    const G4double spin2 = getValue2(y, Value3D::Spin);
    if (spin2 > 0) {
        const G4double spinError2 = getValue2(yError, Value3D::Spin) / spin2;
        error2 = std::max(error2, spinError2);
    }
#endif
    return std::sqrt(error2) / errorTolerance;
}

} // field_utils

