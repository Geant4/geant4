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
// $Id: G4ElectroMagneticField.cc 96678 2016-04-29 16:21:01Z gcosmo $
//
// --------------------------------------------------------------------

#include "G4ElectroMagneticField.hh"

G4ElectroMagneticField::G4ElectroMagneticField()
  : G4Field( false ) // No gravitational field (default)
{
}

G4ElectroMagneticField::~G4ElectroMagneticField()
{
}

G4ElectroMagneticField::G4ElectroMagneticField(const G4ElectroMagneticField &r)
  : G4Field( r )    // To allow extension to joint EM & g field
{
}

G4ElectroMagneticField& 
G4ElectroMagneticField::operator = (const G4ElectroMagneticField &p)
{
  if (&p == this) return *this;
  G4Field::operator=(p);
  return *this;
}
