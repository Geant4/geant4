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
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4VFieldPropagation.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 6 June 2000
// -------------------------------------------------------------------
#include "G4VFieldPropagation.hh"
#include "globals.hh"

G4VFieldPropagation::G4VFieldPropagation()
{ }

G4VFieldPropagation::G4VFieldPropagation(const  G4VFieldPropagation &)
{ }

G4VFieldPropagation::~G4VFieldPropagation()
{ }

const G4VFieldPropagation & G4VFieldPropagation::operator=(const G4VFieldPropagation & )
{
  G4Exception("G4VFieldPropagation::operator= meant not to be accessible");
  return *this;
}

G4int G4VFieldPropagation::operator==(const G4VFieldPropagation & ) const
{
  G4Exception("G4VFieldPropagation::operator== meant not to be accessible");
  return 0;
}

G4int G4VFieldPropagation::operator!=(const G4VFieldPropagation & ) const
{
  G4Exception("G4VFieldPropagation::operator!= meant not to be accessible");
  return 1;
}









