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
// $Id: G4QuadrupoleMagField.cc,v 1.4 2006-06-29 18:24:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "G4QuadrupoleMagField.hh"

G4QuadrupoleMagField::G4QuadrupoleMagField(G4double pGradient)
{
   fGradient = pGradient ;
}

/////////////////////////////////////////////////////////////////////////

G4QuadrupoleMagField::~G4QuadrupoleMagField()
{
}

////////////////////////////////////////////////////////////////////////

void G4QuadrupoleMagField::GetFieldValue( const G4double y[7],
                                                G4double B[3]  ) const  
{
   //   G4double fGradient = 0.001 ;   // Tesla/mm
   B[0] = fGradient*y[1] ;
   B[1] = fGradient*y[0] ;
   B[2] = 0 ;
}
