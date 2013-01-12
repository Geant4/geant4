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
// $Id$
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

#include <fstream>
#include <sstream>

#include "G4HadronicException.hh"
#include "G4SystemOfUnits.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"

using namespace std;

const G4int G4NeutronCaptureXS::amin[] = {0,
					  0, 0, 6, 0,10,12,14,16, 0, 0, //1-10
					  0, 0, 0,28, 0, 0, 0,36, 0,40, //11-20
					  0, 0, 0, 0, 0,54, 0,58,63,64, //21-30
					  0,70, 0, 0, 0, 0, 0, 0, 0,90, //31-40
					  0, 0, 0, 0, 0, 0,107,106, 0,112, //41-50
					  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //51-60
					  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //61-70
					  0, 0, 0,180, 0, 0, 0, 0, 0, 0, //71-80
					  0,204, 0, 0, 0, 0, 0, 0, 0, 0, //81-90
					  0,235};
const G4int G4NeutronCaptureXS::amax[] = {0,
					  0, 0, 7, 0,11,13,15,18, 0, 0, //1-10
					  0, 0, 0,30, 0, 0, 0,40, 0,48, //11-20
					  0, 0, 0, 0, 0,58, 0,64,65,70, //21-30
					  0,76, 0, 0, 0, 0, 0, 0, 0,96, //31-40
					  0, 0, 0, 0, 0, 0,109,116, 0,124, //41-50
					  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //51-60
					  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //61-70
					  0, 0, 0,186, 0, 0, 0, 0, 0, 0, //71-80
					  0,208, 0, 0, 0, 0, 0, 0, 0, 0, //81-90
					  0,238};

G4NeutronCaptureXS::G4NeutronCaptureXS() 
 : G4VCrossSectionDataSet("G4NeutronCaptureXS"),
   emax(20*MeV),maxZ(92)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronCaptureXS::G4NeutronCaptureXS: Initialise for Z < "
	    << maxZ + 1 << G4endl;
  }
  //data.resize(maxZ+1, 0);
  data.SetName("NeutronCapture");
  work.resize(13,0);
  temp.resize(13,0.0);
  isInitialized = false;
}

G4NeutronCaptureXS::~G4NeutronCaptureXS()
{
  /*
  for(G4int i=0; i<=maxZ; ++i) {
    delete data[i];
  }
  */
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
  if(ekin > emax || Z < 1 || Z > maxZ) { return xs; }
  const G4double elimit = 1.0e-10*eV;
  if(ekin < elimit) { ekin = elimit; }

  //  G4PhysicsVector* pv = data[Z];
  G4PhysicsVector* pv = data.GetElementData(Z);

  // element was not initialised
  if(!pv) {
    Initialise(Z);
    //    pv = data[Z];
    pv = data.GetElementData(Z);
    if(!pv) { return xs; }
  }

  G4double e1 = pv->Energy(0);
  if(ekin < e1) { xs = (*pv)[0]*std::sqrt(e1/ekin); }
  else          { xs = pv->Value(ekin); }

  if(verboseLevel > 0){
    G4cout  << "ekin= " << ekin << ",  xs= " << xs << G4endl;
  }
  return xs;
}

G4double 
G4NeutronCaptureXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
				       G4int Z, G4int A,
				       const G4Isotope*, const G4Element*,
				       const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin > emax || Z < 1 || Z > maxZ) { return xs; }
  const G4double elimit = 1.0e-10*eV;
  if(ekin < elimit) { ekin = elimit; }

  //  G4PhysicsVector* pv = data[Z];
  G4PhysicsVector* pv = data.GetElementData(Z);

  // element was not initialised
  if(!pv) {
    Initialise(Z);
    //    pv = data[Z];
    pv = data.GetElementData(Z);
    if(!pv) { return xs; }
  }
  pv = data.GetComponentDataByID(Z, A);
  if(!pv) { return xs; }

  G4double e1 = pv->Energy(0);
  if(ekin < e1) { xs = (*pv)[0]*std::sqrt(e1/ekin); }
  else          { xs = pv->Value(ekin); }

  if(verboseLevel > 0){
    G4cout  << "ekin= " << ekin << ",  xs= " << xs << G4endl;
  }
  return xs;
}

G4Isotope* G4NeutronCaptureXS::SelectIsotope(const G4Element* anElement,
					     G4double kinEnergy)
{
  G4int nIso = anElement->GetNumberOfIsotopes();
  G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
  G4Isotope* iso = (*isoVector)[0];

  // more than 1 isotope
  if(1 < nIso) {
    G4int Z = G4lrint(anElement->GetZ());
    if(Z > maxZ) { Z = maxZ; }
    G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4double q = G4UniformRand();
    G4double sum = 0.0;

    // is there isotope wise cross section?
    if(0 == amin[Z]) {
      for (G4int j = 0; j<nIso; ++j) {
	sum += abundVector[j];
	if(q <= sum) {
	  iso = (*isoVector)[j];
	  break;
	}
      }
    } else {
      size_t nmax = data.GetNumberOfComponents(Z);
      if(temp.size() < nmax) { temp.resize(nmax,0.0); }
      for (size_t i=0; i<nmax; ++i) {
	G4int A = (*isoVector)[i]->GetN();
	G4PhysicsVector* v = data.GetComponentDataByID(Z, A);
        if(v) { sum += abundVector[i]*v->Value(kinEnergy); }
        temp[i] = sum;
      }
      sum *= q;
      for (size_t j = 0; j<nmax; ++j) {
        if(temp[j] >= sum) {
          iso = (*isoVector)[j];
          break;
	}
      }
    }
  }
  return iso;
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
      //      if(!data[Z]) { Initialise(Z, path); }
      if(!data.GetElementData(Z)) { Initialise(Z, path); }
    }
  }
}

void 
G4NeutronCaptureXS::Initialise(G4int Z, const char* p)
{
  if(data.GetElementData(Z)) { return; }
  const char* path = p;

  // check environment variable 
  if(!p) {
    path = getenv("G4NEUTRONXSDATA");
    if (!path) {
      throw G4HadronicException(__FILE__, __LINE__, 
				"G4NEUTRONXSDATA environment variable not defined");
      return;
    }
  }

  // upload element data 
  std::ostringstream ost;
  ost << path << "/cap" << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data.InitialiseForElement(Z, v);

  // upload isotope data
  if(amin[Z] > 0) {
    size_t n = 0;
    size_t i = 0;
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    if(work.size() < nmax) { work.resize(nmax,0); }
    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << path << "/cap" << Z << "_" << A;
      v = RetrieveVector(ost1, false);
      if(v) { ++n; }
      work[i] = v;
      ++i;
    }
    data.InitialiseForComponent(Z, n);
    for(size_t j=0; j<i; ++j) {
      if(work[j]) { data.AddComponent(Z, amin[Z]+j, work[j]); }
    }
  }
}
 
G4PhysicsVector* 
G4NeutronCaptureXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = 0;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(!warn) { return v; }
    G4cout << ost.str() << "  is not opened by G4NeutronCaptureXS" << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,"NO data sets opened");
  }else{
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4NeutronCaptureXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if(!v->Retrieve(filein, true)) {
      G4cout << ost.str() << " is not retrieved in G4NeutronCaptureXS" << G4endl;
      throw G4HadronicException(__FILE__, __LINE__,"ERROR: retrieve data fail");
    }
  }
  return v;
}
