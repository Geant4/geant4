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
// $Id: G4LineCurrentMagField.cc 96751 2016-05-04 09:39:38Z gcosmo $
// -------------------------------------------------------------------

#include "G4LineCurrentMagField.hh"
#include "globals.hh"

G4LineCurrentMagField::G4LineCurrentMagField(G4double pFieldConstant)
{
   fFieldConstant = pFieldConstant ;
}

G4Field* G4LineCurrentMagField::Clone() const
{
    return new G4LineCurrentMagField( fFieldConstant );
}

////////////////////////////////////////////////////////////////////////

G4LineCurrentMagField::~G4LineCurrentMagField()
{
}

///////////////////////////////////////////////////////////////////////////


void G4LineCurrentMagField::GetFieldValue( const G4double yTrack[7],
                                                 G4double B[3]      ) const  
{
   //   G4double fFieldConstant = 100 ;
   G4double a = 1.00 ;   // mm
   G4double x = a*yTrack[0], y = a*yTrack[1] ;
   G4double x2 = x*x, y2 = y*y, r2 = x2 + y2 ;
   G4double r = std::sqrt(r2+a*a) ;
   G4double Br = fFieldConstant/r;
   B[0] = -Br*y/r ;
   B[1] = Br*x/r ;
   B[2] = 0 ;
}

// -----------------------------------------------------------------
