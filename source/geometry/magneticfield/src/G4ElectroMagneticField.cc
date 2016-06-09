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
// $Id: G4ElectroMagneticField.cc,v 1.1 2003/04/02 08:49:11 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// --------------------------------------------------------------------

#include "G4ElectroMagneticField.hh"

G4ElectroMagneticField::G4ElectroMagneticField()
{
}

G4ElectroMagneticField::~G4ElectroMagneticField()
{
}

G4ElectroMagneticField::G4ElectroMagneticField(const G4ElectroMagneticField &r)
  : G4MagneticField(r)
{
}

G4ElectroMagneticField& 
G4ElectroMagneticField::operator = (const G4ElectroMagneticField &p)
{
  if (&p == this) return *this;
  *this = p; return *this;
}
