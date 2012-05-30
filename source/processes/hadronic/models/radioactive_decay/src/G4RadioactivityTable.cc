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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4RadioactivityTable.cc
//
// Version:             0.a
// Date:                29/10/2010
// Author:              F Lei
// Organisation:        QinetiQ UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 October 2010, F. Lei QinetiQ UK
// first release.
//
///////////////////////////////////////////////////////////////////////////////


#include "G4RadioactivityTable.hh"

#include <map>

G4RadioactivityTable::G4RadioactivityTable()
{ 
}
///////////////////////////////////////////////////////////////////////////////
//
G4RadioactivityTable::~G4RadioactivityTable()
{
  fRadioactivity.clear(); 
}
///////////////////////////////////////////////////////////////////////////////
//

G4int G4RadioactivityTable::Entries() const
{
  return (G4int) fRadioactivity.size();
}
///////////////////////////////////////////////////////////////////////////////
//
void G4RadioactivityTable::AddIsotope(G4int Z, G4int A, G4double E, G4double rate, G4double weight)
{
  G4double drate = rate*weight;
  G4double derror = drate*rate;
  G4TwoVector entry = G4TwoVector(drate,derror);
  std::map<G4ThreeVector,G4TwoVector>::iterator it;
  it = fRadioactivity.find(G4ThreeVector(Z,A,E));
  if (it == fRadioactivity.end()) {
    fRadioactivity[G4ThreeVector(Z,A,E)] = entry;
  } else {
    fRadioactivity[G4ThreeVector(Z,A,E)] += entry;
  }
}
///////////////////////////////////////////////////////////////////////////////
//
G4TwoVector G4RadioactivityTable::GetRate(G4int Z, G4int A, G4double E)
{
  if (fRadioactivity.end() == fRadioactivity.find(G4ThreeVector(Z,A,E))) {
    G4cout << G4ThreeVector(Z,A,E) << " is not in the map" << G4endl;
    G4TwoVector rate = G4TwoVector(0.,0.);
    return rate  ;
  }
  else
    return fRadioactivity[G4ThreeVector(Z,A,E)];
}
///////////////////////////////////////////////////////////////////////////////
//
map<G4ThreeVector,G4TwoVector>*  G4RadioactivityTable::GetTheMap()
{
  return &fRadioactivity;
}
///////////////////////////////////////////////////////////////////////////////
//














