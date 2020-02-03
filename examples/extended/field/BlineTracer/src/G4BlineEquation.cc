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
/// \file field/BlineTracer/src/G4BlineEquation.cc
/// \brief Implementation of the G4BlineEquation class
//
//
//
// 
// --------------------------------------------------------------------
//
// G4BlineEquation implementation
//
// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------

#include "G4BlineEquation.hh"

///////////////////////////////////////////////////////////////////////////

G4BlineEquation::G4BlineEquation( G4MagneticField* MagField )
  : G4Mag_EqRhs( MagField ) 
{
  fBackward_direction=false;
  fDirection=1.;
}

///////////////////////////////////////////////////////////////////////////

G4BlineEquation::~G4BlineEquation()
{
}

/////////////////////////////////////////////////////////////////////////////

void G4BlineEquation::EvaluateRhsGivenB( const G4double y[],
                                         const G4double B[3],
                                               G4double dydx[] ) const
{
  G4double Bmag = fDirection*std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  dydx[0] = B[0]/Bmag;       
  dydx[1] = B[1]/Bmag;       
  dydx[2] = B[2]/Bmag;

  dydx[3]=0. * y[0]; //y[0] is used to remove warning
  dydx[4]=0.;
  dydx[5]=0.;
}

//////////////////////////////////////////////////////////////////////

void G4BlineEquation::SetBackwardDirectionOfIntegration(G4bool abool)
{
  fBackward_direction=abool;
  fDirection=1.;
  if (fBackward_direction) fDirection= -1.;
}
