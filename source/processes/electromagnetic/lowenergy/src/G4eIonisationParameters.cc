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
// $Id: G4eIonisationParameters.cc,v 1.2 2001-08-20 17:05:10 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4eIonisationParameters.hh"
#include "G4DataVector.hh"
#include "g4std/fstream"
#include "g4std/strstream"

G4IonisationParameters:: G4IonisationParameters(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ)
{
  LoadData();
}

G4IonisationParameters::~G4IonisationParameters()
{ }

G4double G4IonisationParameters::Parameter(G4int Z, G4int shellIndex, G4int parameterIndex) const
{
  // To be implemented
  G4double value = 0.;
  return value;
}

//const G4DataVector& G4IonisationParameters::Parameters(G4int Z, G4int shellIndex) const
//{
// To be implemented
//}

void G4IonisationParameters::LoadData()
{
  // To be implemented
}

void G4IonisationParameters::PrintData() const
{
  // To be implemented
}
