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
// $Id: G4ElectricField.cc,v 1.1 2003/11/05 10:40:59 japost Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// --------------------------------------------------------------------

#include "G4ElectricField.hh"

G4ElectricField::G4ElectricField()
{
}

G4ElectricField::~G4ElectricField()
{
}

G4ElectricField::G4ElectricField(const G4ElectricField &)
  : G4ElectroMagneticField()
{
}

G4ElectricField& G4ElectricField::operator = (const G4ElectricField &p)
{
  if (&p == this) return *this; *this = p; return *this;
}
