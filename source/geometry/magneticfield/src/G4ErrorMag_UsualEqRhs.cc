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
// $Id: G4ErrorMag_UsualEqRhs.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class implementation file 
// --------------------------------------------------------------------

#include "G4ErrorMag_UsualEqRhs.hh"
#include "G4ErrorPropagatorData.hh"
 
//---------------------------------------------------------------------

G4ErrorMag_UsualEqRhs::G4ErrorMag_UsualEqRhs( G4MagneticField* MagField )
  : G4Mag_UsualEqRhs( MagField )
{
}

G4ErrorMag_UsualEqRhs::~G4ErrorMag_UsualEqRhs()
{
}

//---------------------------------------------------------------------
void
G4ErrorMag_UsualEqRhs::EvaluateRhsGivenB( const G4double y[],
                                          const G4double B[3],
                                                G4double dydx[] ) const
{

  G4Mag_UsualEqRhs::EvaluateRhsGivenB(y, B, dydx );
  
  if(G4ErrorPropagatorData::GetErrorPropagatorData()->GetMode()
     == G4ErrorMode_PropBackwards)
  {
    G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
    G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );
    
    G4double cof = FCof()*inv_momentum_magnitude;

    dydx[3] = cof*(y[4]*(-B[2]) - y[5]*(-B[1])) ;
    dydx[4] = cof*(y[5]*(-B[0]) - y[3]*(-B[2])) ;
    dydx[5] = cof*(y[3]*(-B[1]) - y[4]*(-B[0])) ;    
  }
  return;
}
