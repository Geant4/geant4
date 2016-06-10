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
// $Id: G4NeutronInelasticXS.cc 93682 2015-10-28 10:09:49Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4NeutronInelasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
//

#include "G4NeutronInelasticXS.hh"
#include "G4Neutron.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <sstream>

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4NeutronInelasticXS);

using namespace std;

const G4int G4NeutronInelasticXS::amin[] = {
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
const G4int G4NeutronInelasticXS::amax[] = {
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

G4double G4NeutronInelasticXS::coeff[] = {1.0};

G4ElementData* G4NeutronInelasticXS::data = 0;

G4NeutronInelasticXS::G4NeutronInelasticXS() 
  : G4VCrossSectionDataSet(Default_Name()),
    proton(G4Proton::Proton())
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout << "G4NeutronInelasticXS::G4NeutronInelasticXS Initialise for Z < " 
	    << MAXZINEL << G4endl;
  }
  ggXsection = new G4ComponentGGHadronNucleusXsc();
  fNucleon = new G4HadronNucleonXsc();
  isMaster = false;
}

G4NeutronInelasticXS::~G4NeutronInelasticXS()
{
  //G4cout << "G4NeutronInelasticXS::~G4NeutronInelasticXS() " 
  // << " isMaster= " << isMaster << "  data: " << data << G4endl;
  delete fNucleon;
  delete ggXsection;
  if(isMaster) { delete data; data = 0; }
}

void G4NeutronInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronInelasticXS calculates the neutron inelastic scattering\n"
          << "cross section on nuclei using data from the high precision\n"
          << "neutron database.  These data are simplified and smoothed over\n"
          << "the resonance region in order to reduce CPU time.\n"
          << "G4NeutronInelasticXS is valid for energies up to 20 MeV, for\n"
          << "nuclei through U.\n";
}

G4bool 
G4NeutronInelasticXS::IsElementApplicable(const G4DynamicParticle*, 
					  G4int, const G4Material*)
{
  return true;
}

G4bool 
G4NeutronInelasticXS::IsIsoApplicable(const G4DynamicParticle*,
				      G4int /*ZZ*/, G4int /*AA*/,
				      const G4Element*, const G4Material*)
{
  return true;
}

G4double G4NeutronInelasticXS::GetElementCrossSection(
         const G4DynamicParticle* aParticle,
	 G4int Z, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();

  if(Z >= MAXZINEL) { Z = MAXZINEL - 1; }
  else if(Z < 1)    { Z = 1; }

  G4int Amean = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));

  G4PhysicsVector* pv = data->GetElementData(Z);
  //  G4cout  << "G4NeutronInelasticXS::GetCrossSection e= " << ekin 
  //  << " Z= " << Z << G4endl;

  // element was not initialised
  if(!pv) {
    Initialise(Z);
    pv = data->GetElementData(Z);
    if(!pv) { return xs; }
  }

  G4double e1 = pv->Energy(0);
  if(ekin <= e1) { return xs; }

  G4double e2 = pv->GetMaxEnergy();

  if(ekin <= e2) { 
    xs = pv->Value(ekin); 
  } else if(1 == Z) { 
    fNucleon->GetHadronNucleonXscPDG(aParticle, proton);
    xs = coeff[1]*fNucleon->GetInelasticHadronNucleonXsc();
  } else {          
    ggXsection->GetIsoCrossSection(aParticle, Z, Amean);
    xs = coeff[Z]*ggXsection->GetInelasticGlauberGribovXsc();
  }

  if(verboseLevel > 0) {
    G4cout  << "ekin= " << ekin << ",  XSinel= " << xs << G4endl;
  }
  return xs;
}

G4double G4NeutronInelasticXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int Z, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material*)
{
  if(Z >= MAXZINEL) { Z = MAXZINEL - 1; }
  else if(Z < 1)    { Z = 1; }

  return IsoCrossSection(aParticle->GetKineticEnergy(), Z, A);
}

G4double 
G4NeutronInelasticXS::IsoCrossSection(G4double ekin, G4int Z, G4int A)
{
  G4double xs = 0.0;
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  // element was not initialised
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(!pv) {
    Initialise(Z);
    pv = data->GetElementData(Z);
  }

  // isotope cross section exist
  if(pv && amin[Z] > 0  && A >= amin[Z] && A <= amax[Z]) {
    pv = data->GetComponentDataByIndex(Z, A - amin[Z]);
    if(pv && ekin > pv->Energy(0)) { xs = pv->Value(ekin); }
  }
  if(verboseLevel > 0){
    G4cout  << "ekin= " << ekin << ",  xs= " << xs << G4endl;
  }
  return xs;
}

G4Isotope* G4NeutronInelasticXS::SelectIsotope(
           const G4Element* anElement, G4double kinEnergy)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
  G4Isotope* iso = (*isoVector)[0];

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;

  // more than 1 isotope
  if(1 < nIso) {
    G4int Z = G4lrint(anElement->GetZ());
    //G4cout << "SelectIsotope Z= " << Z << G4endl;

    G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4double q = G4UniformRand();
    G4double sum = 0.0;

    // is there isotope wise cross section?
    size_t j;
    if(0 == amin[Z] || Z >= MAXZINEL) {
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
	//G4cout << j << "-th isotope " << (*isoVector)[j]->GetN() 
	//       <<  " abund= " << abundVector[j] << G4endl;
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
G4NeutronInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0){
    G4cout << "G4NeutronInelasticXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "neutron") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronInelasticXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }

  if(!data) { 
    isMaster = true;
    data = new G4ElementData(); 
    data->SetName("NeutronInelastic");
    temp.resize(13,0.0);
  }

  // it is possible re-initialisation for the new run
  if(isMaster) {

    // check environment variable 
    // Build the complete string identifying the file with the data set
    char* path = getenv("G4NEUTRONXSDATA");

    G4DynamicParticle* dynParticle = 
      new G4DynamicParticle(G4Neutron::Neutron(),G4ThreeVector(1,0,0),1);

    // Access to elements
    const G4ElementTable* theElmTable = G4Element::GetElementTable();
    size_t numOfElm = G4Element::GetNumberOfElements();
    if(numOfElm > 0) {
      for(size_t i=0; i<numOfElm; ++i) {
	G4int Z = G4lrint(((*theElmTable)[i])->GetZ());
	if(Z < 1)              { Z = 1; }
	else if(Z >= MAXZINEL) { Z = MAXZINEL-1; }
	//G4cout << "Z= " << Z << G4endl;
	// Initialisation 
	if(!(data->GetElementData(Z))) { 
	  Initialise(Z, dynParticle, path); 
	}
      }
    }
    delete dynParticle;
  }
}

void 
G4NeutronInelasticXS::Initialise(G4int Z, G4DynamicParticle* dp, 
				 const char* p)
{
  if(data->GetElementData(Z) || Z < 1 || Z >= MAXZINEL) { return; }
  const char* path = p;
  if(!p) {
    // check environment variable 
    // Build the complete string identifying the file with the data set
    path = getenv("G4NEUTRONXSDATA");
    if (!path) {
      G4Exception("G4NeutronInelasticXS::Initialise(..)","had013",
		  FatalException,
                  "Environment variable G4NEUTRONXSDATA is not defined");
      return;
    }
  }
  G4DynamicParticle* dynParticle = dp;
  if(!dp) {
    dynParticle = 
      new G4DynamicParticle(G4Neutron::Neutron(),G4ThreeVector(1,0,0),1);
  }

  G4int Amean = G4lrint(G4NistManager::Instance()->GetAtomicMassAmu(Z));

  // upload element data 
  std::ostringstream ost;
  ost << path << "/inelast" << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);
  /*
  G4cout << "G4NeutronInelasticXS::Initialise for Z= " << Z 
	 << " A= " << Amean << "  Amin= " << amin[Z] 
	 << "  Amax= " << amax[Z] << G4endl;
  */
  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << path << "/inelast" << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v1); 
    }
  }

  // smooth transition 
  G4double emax = v->GetMaxEnergy();
  G4double sig1 = (*v)[v->GetVectorLength() - 1];
  dynParticle->SetKineticEnergy(emax);
  G4double sig2 = 0.0;
  if(1 == Z) {
    fNucleon->GetHadronNucleonXscPDG(dynParticle, proton);
    sig2 = fNucleon->GetInelasticHadronNucleonXsc();
  } else {
    ggXsection->GetIsoCrossSection(dynParticle, Z, Amean);
    sig2 = ggXsection->GetInelasticGlauberGribovXsc();
  }
  if(sig2 > 0.) { coeff[Z] = sig1/sig2; } 
  if(!dp) { delete dynParticle; }
}

G4PhysicsVector* 
G4NeutronInelasticXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = 0;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(warn) { 
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4NEUTRONXSDATA");
    }
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4NeutronInelasticXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4NEUTRONXSDATA");
    }
  }
  return v;
}
