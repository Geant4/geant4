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
//      File name:     G4VNuclearField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------
#include "G4VNuclearField.hh"
#include "globals.hh"

G4VNuclearField::G4VNuclearField(G4V3DNucleus * aNucleus) : 
theNucleus(aNucleus),
radius(aNucleus->GetOuterRadius() + 4*fermi)
{
}

G4VNuclearField::G4VNuclearField(const  G4VNuclearField &right) :
theNucleus(right.theNucleus),
radius(right.radius)
{
}

G4VNuclearField::~G4VNuclearField()
{
}

const G4VNuclearField & G4VNuclearField::operator=(const G4VNuclearField & )
{
  G4Exception("G4VNuclearField::operator= meant not to be accessible");
  return *this;
}

G4int G4VNuclearField::operator==(const G4VNuclearField & ) const
{
  G4Exception("G4VNuclearField::operator== meant not to be accessible");
  return 0;
}

G4int G4VNuclearField::operator!=(const G4VNuclearField & ) const
{
  G4Exception("G4VNuclearField::operator!= meant not to be accessible");
  return 1;
}









