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
#include "G4FieldPropagation.hh"
#include "G4HadronicException.hh"

const G4FieldPropagation & G4FieldPropagation::operator=(const G4FieldPropagation &)
{
   throw G4HadronicException(__FILE__, __LINE__, "G4FieldPropagation::operator= meant to be private");
   return *this;
}

int G4FieldPropagation::operator==(const G4FieldPropagation &) const
{
   return 1;
}

int G4FieldPropagation::operator!=(const G4FieldPropagation &) const
{
   return 0;
}
