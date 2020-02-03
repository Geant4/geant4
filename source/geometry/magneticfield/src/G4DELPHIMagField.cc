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
// G4DELPHIMagField implementation
//
// Created: V.Grichine - 03.02.1997
// -------------------------------------------------------------------

#include "G4DELPHIMagField.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

G4DELPHIMagField::G4DELPHIMagField()
{
}

////////////////////////////////////////////////////////////////////////

G4DELPHIMagField::~G4DELPHIMagField()
{
}

///////////////////////////////////////////////////////////////////////

G4Field* G4DELPHIMagField::Clone() const
{
   return new G4DELPHIMagField;
}

///////////////////////////////////////////////////////////////////////

void G4DELPHIMagField::GetFieldValue( const G4double yTrack[7],
                                            G4double B[3]     ) const 
{
   G4int i, n = 8 ;
   G4double a = 0.001 ;   // mm -> m 
   G4double x = a*yTrack[0], y = a*yTrack[1], z = a*yTrack[2] ;
   G4double x2 = x*x, y2 = y*y, z2 = z*z, r2 = x2 + y2 ;
   G4double r4 = r2*r2, z4 = z2*z2, r6 = r4*r2, z6 = z4*z2 ;
   G4double r8 = r4*r4, z8 = z4*z4, r10 = r8*r2, z10 = z8*z2 ;
   G4double rz = z*std::sqrt(r2), r = std::sqrt(r2+a*a) ;
   G4double Br ;
   G4double P[8], Q[8] ; 
   static G4ThreadLocal G4double c[8] =
                          {
                            -9.26e-5, -3.51e-5, 2.94e-6, -1.10e-6, 
                             6.25e-8, -1.77e-8, -6.88e-10, -7.52e-11 
                          } ;
   P[0] = 2*rz ;
   P[1] = 4*rz*(r2 - 4.0*z2/3.0) ;
   P[2] = 6*rz*(r4 - 4*r2*z2 + 1.6*z4) ;
   P[3] = 8*rz*(r6 - 8*r4*z2 + 9.6*r2*z4 - 64.0*z6/35.0) ;
   P[4] = 10*rz*(r8 - 40.0*r6*z2/3.0 + 32*r4*z4
                    - 128.0*r2*z6/7.0 + 128.0*z8/63.0);
   P[5] = 0 ;
   P[6] = 0 ;
   P[7] = 0 ;
   
   Q[0] = r2 - 2*z2 ;
   Q[1] = r4 - 8*r2*z2 + 8.0*z4/3.0 ;
   Q[2] = r6 - 18*r4*z2 + 24*r2*z4 - 3.2*z6 ;
   Q[3] = r8 - 32*r6*z2 + 96*r4*z4 - 51.2*r2*z6 +128.0*z8/35.0 ;
   Q[4] = r10 - 50*r8*z2 + 800.0*r6*z4/3.0 - 320*r4*z6 
          + 640.0*r2*z8/7.0 - 256.0*z10/63.0 ;
   Q[5] = 0 ;
   Q[6] = 0 ;
   Q[7] = 0 ;
   
   Br = 0 ;
   B[2] = 1.2*tesla  ;    //  the principal Bz value of DELPHI detector
   for(i=0; i<n; ++i)
   {
      Br += c[i]*P[i] ;
      B[2] += c[i]*Q[i] ;
   }
   B[0] = Br*x/r ;
   B[1] = Br*y/r ;
   return ;
}

// -----------------------------------------------------------------
