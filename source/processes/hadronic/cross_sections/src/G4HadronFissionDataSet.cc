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
// $Id: G4HadronFissionDataSet.cc 68720 2013-04-05 09:18:58Z gcosmo $
//
//
// G4 Physics class: HadronFissionDataSet for cross sections
// F.W. Jones, TRIUMF, 19-MAY-98

#include "G4HadronFissionDataSet.hh"
#include "G4DynamicParticle.hh"
//#include "G4NistManager.hh"
//#include "G4HadTmpUtil.hh"
#include <iostream>


G4HadronFissionDataSet::G4HadronFissionDataSet(const G4String& nam)
  : G4VCrossSectionDataSet(nam)
{
  //theHadronCrossSections = G4HadronCrossSections::Instance();
}


G4HadronFissionDataSet::~G4HadronFissionDataSet() {}


void G4HadronFissionDataSet::CrossSectionDescription(std::ostream& outFile) const 
{
  outFile << "G4HadronFissionDataSet contains cross sections for\n"
          << "neutron-induced fission of nuclei.  They were developed as\n"
          << "part of the Gheisha hadronic package by H. Fesefeldt.  The\n"
          << "cross sections are valid for all incident neutron energies.\n";
}


G4bool
G4HadronFissionDataSet::IsElementApplicable(const G4DynamicParticle*, 
					    G4int,
					    const G4Material*)
{
  return true;
  //  return theHadronCrossSections->IsApplicable(aParticle);
}

G4double
G4HadronFissionDataSet::GetElementCrossSection(const G4DynamicParticle*, 
					       G4int, 
					       const G4Material*)
{
  return 0.0;
}

/*
G4double
G4HadronFissionDataSet::GetElementCrossSection(const G4DynamicParticle* aParticle, 
					       G4int Z, 
					       const G4Material*)
{
  G4int A = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));
  return theHadronCrossSections->GetFissionCrossSection(aParticle, Z, A);
}
*/
