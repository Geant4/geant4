//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4BlineEquation.cc,v 1.1 2003/11/25 09:29:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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
  backward_direction=false;
  direction=1.;
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
  G4double Bmag = direction*std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
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
  backward_direction=abool;
  direction=1.;
  if (backward_direction) direction= -1.;
}
