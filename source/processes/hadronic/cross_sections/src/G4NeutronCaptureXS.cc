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
// $Id: G4NeutronCaptureXS.cc,v 1.8 2011-01-09 02:37:48 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4NeutronCaptureXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
//

#include "G4HadronicException.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4DynamicParticle.hh"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

G4NeutronCaptureXS::G4NeutronCaptureXS() 
 : G4VCrossSectionDataSet("G4NeutronCaptureXS"),
   emax(20*MeV),maxZ(92)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronCaptureXS::G4NeutronCaptureXS: Initialise for Z < "
	    << maxZ + 1 << G4endl;
  }
  data.resize(maxZ+1, 0);
  isInitialized = false;
}

G4NeutronCaptureXS::~G4NeutronCaptureXS()
{
  for(G4int i=0; i<=maxZ; ++i) {
    delete data[i];
  }
}

void G4NeutronCaptureXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronCaptureXS calculates the neutron capture cross sections\n"
          << "on nuclei using data from the high precision neutron database.\n"
          << "These data are simplified and smoothed over the resonance region\n"
          << "in order to reduce CPU time.  G4NeutronCaptureXS is valid up to\n"
          << "20 MeV for all targets through U.\n";
}
 
G4bool 
G4NeutronCaptureXS::IsElementApplicable(const G4DynamicParticle*, 
					G4int, const G4Material*)
{
  return true;
}

G4bool 
G4NeutronCaptureXS::IsIsoApplicable(const G4DynamicParticle*,
				    G4int /*ZZ*/, G4int /*AA*/,
				    const G4Element*, const G4Material*)
{
  return false;
}

G4double 
G4NeutronCaptureXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
					   G4int Z, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin > emax) { return xs; }
  const G4double elimit = 1.0e-10*eV;
  if(ekin < elimit) { ekin = elimit; }

  if(Z < 1 || Z > maxZ) { return xs; }
  G4PhysicsVector* pv = data[Z];

  // element was not initialised
  if(!pv) {
    Initialise(Z);
    pv = data[Z];
    if(!pv) { return xs; }
  }

  G4int n = pv->GetVectorLength() - 1;
  G4double e1 = pv->Energy(0);
  G4double e2 = pv->Energy(n);
  if(ekin < e1)       { xs = (*pv)[0]*std::sqrt(e1/ekin); }
  else if(ekin <= e2) { xs = pv->Value(ekin); }

  if(verboseLevel > 0){
    G4cout  << "ekin= " << ekin << ",  xs= " << xs << G4endl;
  }
  return xs;
}

void 
G4NeutronCaptureXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(isInitialized) { return; }
  if(verboseLevel > 0){
    G4cout << "G4NeutronCaptureXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "neutron") { 
    throw G4HadronicException(__FILE__, __LINE__,"Wrong particle type");
    return; 
  }
  isInitialized = true;

  // check environment variable 
  // Build the complete string identifying the file with the data set
  char* path = getenv("G4NEUTRONXSDATA");
  if (!path){
    throw G4HadronicException(__FILE__, __LINE__, 
			      "G4NEUTRONXSDATA environment variable not defined");
    return;
  }

  // Access to elements
  const G4ElementTable* theElmTable = G4Element::GetElementTable();
  size_t numOfElm = G4Element::GetNumberOfElements();
  if(numOfElm > 0) {
    for(size_t i=0; i<numOfElm; ++i) {
      G4int Z = G4int(((*theElmTable)[i])->GetZ());
      if(Z < 1)         { Z = 1; }
      else if(Z > maxZ) { Z = maxZ; }
      //G4cout << "Z= " << Z << G4endl;
      // Initialisation 
      if(!data[Z]) { Initialise(Z, path); }
    }
  }
}

void 
G4NeutronCaptureXS::Initialise(G4int Z, const char* p)
{
  if(data[Z]) { return; }
  const char* path = p;
  if(!p) {
  // check environment variable 
  // Build the complete string identifying the file with the data set
    path = getenv("G4NEUTRONXSDATA");
    if (!path) {
      throw G4HadronicException(__FILE__, __LINE__, 
				"G4NEUTRONXSDATA environment variable not defined");
      return;
    }
  }

  // upload data from file
  data[Z] = new G4PhysicsLogVector();

  std::ostringstream ost;
  ost << path << "/cap" << Z ;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    G4cout << ost.str() << "  is not opened by G4NeutronCaptureXS" << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,"NO data sets opened");
    return;
  }else{
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4NeutronCaptureXS" << G4endl;
    }
    // retrieve data from DB
    data[Z]->Retrieve(filein, true);
  } 
}
 
