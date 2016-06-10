//
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
// $Id: G4ErrorMag_UsualEqRhs.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class header file 
// --------------------------------------------------------------------
//
// Class description:
//
// Serves to reverse the magnetic field when propagation is backwards
// for error propagation.

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorMag_UsualEqRhs_hh
#define G4ErrorMag_UsualEqRhs_hh

#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"

class G4ErrorMag_UsualEqRhs : public G4Mag_UsualEqRhs
{
   public:  // with description

     G4ErrorMag_UsualEqRhs( G4MagneticField* MagField );
    ~G4ErrorMag_UsualEqRhs();

     void EvaluateRhsGivenB( const G4double y[],
                             const G4double B[3],
                                   G4double dydx[] ) const;
       // Reverses dedx if propagation is backwards
};

#endif /* G4MAG_USUAL_EQRHS */
