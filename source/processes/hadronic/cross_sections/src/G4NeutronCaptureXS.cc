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
// $Id: G4NeutronCaptureXS.cc 91903 2015-08-10 12:10:36Z gcosmo $
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

#include "G4SystemOfUnits.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4NeutronCaptureXS);

using namespace std;

const G4int G4NeutronCaptureXS::amin[] = {
  0,
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
const G4int G4NeutronCaptureXS::amax[] = {
  0,
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

G4ElementData* G4NeutronCaptureXS::data = 0;

G4NeutronCaptureXS::G4NeutronCaptureXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   emax(20*MeV),elimit(1.0e-10*eV)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronCaptureXS::G4NeutronCaptureXS: Initialise for Z < "
	    << MAXZCAPTURE << G4endl;
  }
  isMaster = false;
}

G4NeutronCaptureXS::~G4NeutronCaptureXS()
{
  if(isMaster) { delete data; data = 0; }
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
  return true;
}

G4double 
G4NeutronCaptureXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
					   G4int Z, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin > emax || Z < 1 || Z >= MAXZCAPTURE) { return xs; }
  if(ekin < elimit) { ekin = elimit; }

  // element was not initialised
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(!pv) {
    Initialise(Z);
    pv = data->GetElementData(Z);
    if(!pv) { return xs; }
  }

  G4double e1 = pv->Energy(0);
  if(ekin < e1) { xs = (*pv)[0]*std::sqrt(e1/ekin); }
  else if(ekin <= pv->GetMaxEnergy()) { xs = pv->Value(ekin); }

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
  if(ekin <= emax && Z > 0 && Z < MAXZCAPTURE) { 
    xs = IsoCrossSection(ekin, Z, A); 
  }
  return xs;
}

G4double G4NeutronCaptureXS::IsoCrossSection(G4double ekin, G4int Z, G4int A)
{
  G4double xs = 0.0;
  if(ekin < elimit) { ekin = elimit; }

  // element was not initialised
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(!pv) {
    Initialise(Z);
    pv = data->GetElementData(Z);
  }

  // isotope cross section exist
  if(pv && amin[Z] > 0 && A >= amin[Z] && A <= amax[Z]) {
    pv  = data->GetComponentDataByID(Z, A - amin[Z]);
    if(pv) { 

      G4double e1 = pv->Energy(1);
      if(ekin < e1) { xs = (*pv)[1]*std::sqrt(e1/ekin); }
      else if(ekin <= pv->GetMaxEnergy()) { xs = pv->Value(ekin); }
    }
  }
  if(verboseLevel > 0) {
    G4cout  << "G4NeutronCaptureXS::IsoCrossSection: Ekin(MeV)= " << ekin/MeV 
	    << "  xs(b)= " << xs/barn 
	    << "  Z= " << Z << "  A= " << A << G4endl;
  }
  return xs;
}

G4Isotope* G4NeutronCaptureXS::SelectIsotope(const G4Element* anElement,
					     G4double kinEnergy)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
  G4Isotope* iso = (*isoVector)[0];

  // more than 1 isotope
  if(1 < nIso) {
    G4int Z = G4lrint(anElement->GetZ());

    G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4double q = G4UniformRand();
    G4double sum = 0.0;

    // is there isotope wise cross section?
    size_t j;
    if(0 == amin[Z] || Z >= MAXZCAPTURE) {
      for (j = 0; j<nIso; ++j) {
	sum += abundVector[j];
	if(q <= sum) {
	  iso = (*isoVector)[j];
	  break;
	}
      }
    } else {

      // element may be not initialised in unit test
      if(!data->GetElementData(Z)) { Initialise(Z); }
      size_t nn = temp.size();
      if(nn < nIso) { temp.resize(nIso, 0.); }

      for (j=0; j<nIso; ++j) {
        sum += abundVector[j]*IsoCrossSection(kinEnergy, Z, 
					      (*isoVector)[j]->GetN());
        temp[j] = sum;
      }
      sum *= q;
      for (j = 0; j<nIso; ++j) {
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
  if(verboseLevel > 0){
    G4cout << "G4NeutronCaptureXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "neutron") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronCaptureXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }

  if(!data) { 
    isMaster = true;
    data = new G4ElementData(); 
    data->SetName("NeutronCapture");
    temp.resize(13,0.0);
  }

  // it is possible re-initialisation for the second run
  if(isMaster) {

    // check environment variable 
    // Build the complete string identifying the file with the data set
    char* path = getenv("G4NEUTRONXSDATA");

    // Access to elements
    const G4ElementTable* theElmTable = G4Element::GetElementTable();
    size_t numOfElm = G4Element::GetNumberOfElements();
    if(numOfElm > 0) {
      for(size_t i=0; i<numOfElm; ++i) {
	G4int Z = G4int(((*theElmTable)[i])->GetZ());
	if(Z < 1)                 { Z = 1; }
	else if(Z >= MAXZCAPTURE) { Z = MAXZCAPTURE-1; }
	//G4cout << "Z= " << Z << G4endl;
	// Initialisation 
	if(!data->GetElementData(Z)) { Initialise(Z, path); }
      }
    }
  }
}

void 
G4NeutronCaptureXS::Initialise(G4int Z, const char* p)
{
  if(data->GetElementData(Z) || Z < 1 || Z >= MAXZCAPTURE) { return; }
  const char* path = p;

  // check environment variable 
  if(!p) {
    path = getenv("G4NEUTRONXSDATA");
    if (!path) {
      G4Exception("G4NeutronCaptureXS::Initialise(..)","had013",FatalException,
                  "Environment variable G4NEUTRONXSDATA is not defined");
      return;
    }
  }

  // upload element data 
  std::ostringstream ost;
  ost << path << "/cap" << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);

  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << path << "/cap" << Z << "_" << A;
      v = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v); 
    }
  }
}
 
G4PhysicsVector* 
G4NeutronCaptureXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = 0;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(warn) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4NeutronCaptureXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4NEUTRONXSDATA");
    }
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4NeutronCaptureXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4NeutronCaptureXS::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4NEUTRONXSDATA");
    }
  }
  return v;
}
