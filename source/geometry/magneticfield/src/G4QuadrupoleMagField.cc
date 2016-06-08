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
#include "G4QuadrupoleMagField.hh"
#include "globals.hh"
#include "geomdefs.hh"



G4QuadrupoleMagField::G4QuadrupoleMagField(G4double pGradient)
{
   fGradient = pGradient ;
}

/////////////////////////////////////////////////////////////////////////

G4QuadrupoleMagField::~G4QuadrupoleMagField()
{
   ;
}

////////////////////////////////////////////////////////////////////////


void
   G4QuadrupoleMagField::GetFieldValue( const G4double y[7],
                                              G4double B[3]  ) const  
{
   //   G4double fGradient = 0.001 ;   // Tesla/mm
   B[0] = fGradient*y[1] ;
   B[1] = fGradient*y[0] ;
   B[2] = 0 ;
   return ;
}



