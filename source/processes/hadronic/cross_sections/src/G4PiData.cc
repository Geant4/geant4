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
// --------------------------------------------------------------------
// by J.P Wellisch, Sun Sep 15 2002.

#include "G4PiData.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"

///////////////////////////////////////////////////////////////////////

G4PiData::G4PiData(const G4double * aT, const G4double * aIn, 
                   const G4double * anE, G4int nP)
{
  G4int i=0;

  for( i = 0; i < nP; i++ )
  {
    std::pair<G4double, G4double> x;
    x.first=aT[i]*millibarn;
    x.second=aIn[i]*millibarn;
    std::pair<G4double, std::pair<G4double, G4double > > aP;
    aP.first=anE[i]*GeV;
    aP.second=x;
    push_back(aP);
  }
}

////////////////////////////////////////////////////////////////////////

G4bool G4PiData::AppliesTo(G4double kineticEnergy)
{
  G4bool result = true;
  if(kineticEnergy>back().first) result = false;
  return result;
}

//////////////////////////////////////////////////////////////////////////

G4double G4PiData::ReactionXSection(G4double kineticEnergy)
{
  G4double result = 0;
  G4PiData::iterator it=begin();
  while(it!=end()&&kineticEnergy>(*it).first) {it++;}  /* Loop checking, 08.01.2016, W. Pokorski */
  if(it==end()) 
  {
    throw G4HadronicException(__FILE__, __LINE__,
        "G4PiData::ReactionXSection: used outside validity range");
  }
  if(it==begin()) it++;
  G4double x1,x2,e1,e2;
  e1=(*(it-1)).first;
  x1=(*(it-1)).second.second;
  e2=(*(it)).first;
  x2=(*(it)).second.second;
  result = std::max(0., x1 + (kineticEnergy-e1)*(x2-x1)/(e2-e1));
  return result;
}

////////////////////////////////////////////////////////////////////////////

G4double G4PiData::ElasticXSection(G4double kineticEnergy)
{
  G4double result = 0;
  G4PiData::iterator it=begin();
  while(it!=end()&&kineticEnergy>(*it).first) {it++;}  /* Loop checking, 08.01.2016, W. Pokorski */
  if(it==end()) 
  {
    throw G4HadronicException(__FILE__, __LINE__,
        "G4PiData::ElasticXSection: used outside validity range");
  }
  if(it==begin()) it++;
  G4double x1,x2,e1,e2;
  e1=(*(it-1)).first;
  x1=(*(it-1)).second.first - (*(it-1)).second.second;
  e2=(*(it)).first;
  x2=(*(it)).second.first - (*(it)).second.second;
  result = std::max(0., x1 + (kineticEnergy-e1)*(x2-x1)/(e2-e1));
  return result;
}

////////////////////////////////////////////////////////////////////////////

G4double G4PiData::TotalXSection(G4double kineticEnergy)
{
  G4double result = 0;
  G4PiData::iterator it=begin();
  while(it!=end()&&kineticEnergy>(*it).first) {it++;}  /* Loop checking, 08.01.2016, W. Pokorski */
  if(it==end()) 
  {
    throw G4HadronicException(__FILE__, __LINE__,
        "G4PiData::TotalXSection: used outside validity range");
  }
  if(it==begin()) it++;
  G4double x1,x2,e1,e2;
  e1=(*(it-1)).first;
  x1=(*(it-1)).second.first;
  e2=(*(it)).first;
  x2=(*(it)).second.first;
  result = std::max(0., x1 + (kineticEnergy-e1)*(x2-x1)/(e2-e1));
  return result;
}
