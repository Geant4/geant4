/
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
// $Id: G4LineCurrentMagField.cc,v 1.1 2002-03-28 13:36:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4LineCurrentMagField.hh"
#include "globals.hh"
#include "geomdefs.hh"



G4LineCurrentMagField::G4LineCurrentMagField(G4double pFieldConstant)
{
   fFieldConstant = pFieldConstant ;
}
// --------------------------------------------------------------------------

G4LineCurrentMagField::~G4LineCurrentMagField()
{
   ;
}

// ------------------------------------------------------------------------


void
   G4LineCurrentMagField::MagneticField (const G4double yTrack[7],
                                         G4double B[3]         )  
{
   //   G4double fFieldConstant = 100 ;
   G4double a = 1.00 ;   // mm -> m 
   G4double x = a*yTrack[0], y = a*yTrack[1], z = a*yTrack[2] ;
   G4double x2 = x*x, y2 = y*y, r2 = x2 + y2 ;
   G4double r = sqrt(r2+a*a) ;
   G4double Br = fFieldConstant/r;
   B[0] = -Br*y/r ;
   B[1] = Br*x/r ;
   B[2] = 0 ;
   return ;
}

// -----------------------------------------------------------------

