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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4ParticleInelasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
//

#include "G4ParticleInelasticXS.hh"
#include "G4Neutron.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4HadronNucleonXsc.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

const G4int G4ParticleInelasticXS::amin[] = {
  0,
  1,   4,   6,   9, 10,  12,  14,  16,  19,  20,  //1-10
 23,  24,  27,  28, 31,  32,  35,  36,  39,  40,  //11-20
 45,  46,  50,  50, 55,  54,  59,  58,  63,  64,  //21-30
 69,  70,  75,   0,  0,   0,   0,   0,   0,  90,  //31-40
  0,  92,   0,   0,  0, 102, 107, 106, 113, 112,  //41-50
  0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  //51-60
  0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  //61-70
  0,   0, 181, 180,  0,   0,   0, 192, 197,   0,  //71-80
  0, 204, 209,   0,  0,   0,   0,   0,   0,   0,  //81-90
  0, 235};
const G4int G4ParticleInelasticXS::amax[] = {
  0,
  1,   4,   7,   9, 11,  13,  15,  18,  19,  22,  //1-10
 23,  26,  27,  30, 31,  34,  37,  40,  41,  48,  //11-20
 45,  50,  51,  54, 55,  58,  59,  64,  65,  70,  //21-30
 71,  76,  75,   0,  0,   0,   0,   0,   0,  96,  //31-40
  0, 100,   0,   0,  0, 110, 109, 116, 115, 124,  //41-50
  0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  //51-60
  0,   0,   0,   0,  0,   0,   0,   0,   0,   0,  //61-70
  0,   0, 181, 186,  0,   0,   0, 198, 197,   0,  //71-80
  0, 208, 209,   0,  0,   0,   0,   0,   0,   0,  //81-90
  0, 238};

G4double G4ParticleInelasticXS::coeff[] = {1.0};

G4ElementData* G4ParticleInelasticXS::data = nullptr;

#ifdef G4MULTITHREADED
  G4Mutex G4ParticleInelasticXS::particleInelasticXSMutex = G4MUTEX_INITIALIZER;
#endif

G4ParticleInelasticXS::G4ParticleInelasticXS(const G4ParticleDefinition* part) 
  : G4VCrossSectionDataSet("G4ParticleInelasticXS"),
    ggXsection(nullptr),
    nnXsection(nullptr),
    fNucleon(nullptr),
    particle(part),
    proton(G4Proton::Proton()),
    particleName(""),
    isMaster(false),
    emax(20*CLHEP::MeV)
{
  if(!part) {
    G4Exception("G4ParticleInelasticXS::G4ParticleInelasticXS(..)","had015",
		FatalException, "NO particle definition in constructor");
  } else {
    //  verboseLevel = 0;
    particleName = particle->GetParticleName();
    if(verboseLevel > 0){
      G4cout << "G4ParticleInelasticXS::G4ParticleInelasticXS for " 
	     << particleName << " on atoms with Z < " << MAXZINELP << G4endl;
    }
    if(particleName == "neutron" || particleName == "proton") {
      ggXsection = new G4ComponentGGHadronNucleusXsc();
      fNucleon = new G4HadronNucleonXsc();
    } else {
      nnXsection = new G4ComponentGGNuclNuclXsc();
    }
  }
  SetForAllAtomsAndEnergies(true);
  fNist  = G4NistManager::Instance();
}

G4ParticleInelasticXS::~G4ParticleInelasticXS()
{
  //G4cout << "G4ParticleInelasticXS::~G4ParticleInelasticXS() " 
  // << " isMaster= " << isMaster << "  data: " << data << G4endl;
  delete fNucleon;
  if(isMaster) { delete data; data = nullptr; }
}

void G4ParticleInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ParticleInelasticXS calculates n, p, d, t, he3, he4 inelastic \n"
          << "cross section on nuclei using data from the high precision\n"
          << "neutron database.  These data are simplified and smoothed over\n"
          << "the resonance region in order to reduce CPU time.\n"
          << "For high energy Glauber-Gribov cross section model is used\n";
}

G4bool 
G4ParticleInelasticXS::IsElementApplicable(const G4DynamicParticle*, 
					   G4int, const G4Material*)
{
  return true;
}

G4bool 
G4ParticleInelasticXS::IsIsoApplicable(const G4DynamicParticle*,
				       G4int, G4int,
				       const G4Element*, const G4Material*)
{
  return true;
}

G4double G4ParticleInelasticXS::GetElementCrossSection(
         const G4DynamicParticle* aParticle,
	 G4int ZZ, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();

  G4int Z = (ZZ >= MAXZINELP) ? MAXZINELP - 1 : ZZ; 

  G4PhysicsVector* pv = data->GetElementData(Z);
  //  G4cout  << "G4ParticleInelasticXS::GetCrossSection e= " << ekin 
  //  << " Z= " << Z << G4endl;

  // element was not initialised or below threshold
  if(!pv || ekin <= pv->Energy(0)) { return xs; }

  if(ekin <= pv->GetMaxEnergy()) { 
    xs = pv->Value(ekin); 
  } else if(1 == Z) {
    if(fNucleon) { 
      fNucleon->GetHadronNucleonXscNS(aParticle, proton);
      xs = coeff[1]*fNucleon->GetInelasticHadronNucleonXsc();
    } else {
      nnXsection->GetZandACrossSection(aParticle, 1, 1);
      xs = coeff[1]*nnXsection->GetInelasticGlauberGribovXsc();
    }
  } else {
    G4int Amean = G4lrint(fNist->GetAtomicMassAmu(Z));
    if(ggXsection) {
      ggXsection->GetIsoCrossSection(aParticle, Z, Amean);
      xs = coeff[Z]*ggXsection->GetInelasticGlauberGribovXsc();
    } else {
      nnXsection->GetZandACrossSection(aParticle, Z, Amean);
      xs = coeff[Z]*nnXsection->GetInelasticGlauberGribovXsc();
    }
  }

  if(verboseLevel > 1) {
    G4cout  << "ElmXS: Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << " xs(bn)= " << xs/CLHEP::barn << " element data for "
	    << particleName << G4endl;
  }
  return xs;
}

G4double G4ParticleInelasticXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int Z, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), Z, A);
}

G4double 
G4ParticleInelasticXS::IsoCrossSection(G4double ekin, G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  G4int Z = (ZZ >= MAXZINELP) ? MAXZINELP - 1 : ZZ; 
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  // first compute isotope cross section 
  if(ekin <=emax && amin[Z]>0 && A >= amin[Z] && A <= amax[Z]) {
    G4PhysicsVector* pviso = data->GetComponentDataByIndex(Z, A - amin[Z]);
    if(pviso) { 
      xs = pviso->Value(ekin); 
      if(verboseLevel > 0){
	G4cout  << "IsoXS for " << particleName 
		<< " Target Z= " << Z << " A= " << A 
		<< " Ekin(MeV)= " << ekin/CLHEP::MeV 
		<< " xs(bn)= " << xs/CLHEP::barn << G4endl;
      }
      return xs;
    }
  }
  // isotope data are not available or applicable
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(pv) { xs = pv->Value(ekin); }
  if(verboseLevel > 0) {
    G4cout  << "IsoXS for " << particleName 
	    << " Target Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << " xs(bn)= " << xs/CLHEP::barn << G4endl;
  }
  return xs;
}

const G4Isotope* G4ParticleInelasticXS::SelectIsotope(
           const G4Element* anElement, G4double kinEnergy)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  // more than 1 isotope
  G4int Z = anElement->GetZasInt();
  //G4cout << "SelectIsotope Z= " << Z << G4endl;

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;

  // is there isotope wise cross section?
  size_t j;
  if(kinEnergy > emax || 0 == amin[Z] || Z >= MAXZINELP) {
    for (j = 0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = anElement->GetIsotope(j);
	break;
      }
    }
    return iso;
  }

  size_t nn = temp.size();
  if(nn < nIso) { temp.resize(nIso, 0.); }

  for (j=0; j<nIso; ++j) {
    //G4cout << j << "-th isotope " << (*isoVector)[j]->GetN() 
    //       <<  " abund= " << abundVector[j] << G4endl;
    sum += abundVector[j]*IsoCrossSection(kinEnergy, Z, 
					  anElement->GetIsotope(j)->GetN());
    temp[j] = sum;
  }
  sum *= q;
  for (j = 0; j<nIso; ++j) {
    if(temp[j] >= sum) {
      iso = anElement->GetIsotope(j);
      break;
    }
  }
  return iso;
}

void 
G4ParticleInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0){
    G4cout << "G4ParticleInelasticXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(&p != particle) { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << particleName << " is expected";
    G4Exception("G4ParticleInelasticXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }

  if(!data) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&particleInelasticXSMutex);
    if(!data) { 
#endif
      isMaster = true;
      data = new G4ElementData(); 
      data->SetName(particleName + "Inelastic");
      temp.resize(13,0.0);
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&particleInelasticXSMutex);
#endif
  }

  // it is possible re-initialisation for the new run
  if(isMaster) {

    // check environment variable 
    // Build the complete string identifying the file with the data set
    char* path = getenv("G4PARTICLEXSDATA");

    G4DynamicParticle* dynParticle = 
      new G4DynamicParticle(particle,G4ThreeVector(1,0,0),1);

    // Access to elements
    const G4ElementTable* theElmTable = G4Element::GetElementTable();
    size_t numOfElm = G4Element::GetNumberOfElements();
    for(size_t i=0; i<numOfElm; ++i) {
      G4int Z = ((*theElmTable)[i])->GetZasInt();
      if(Z >= MAXZINELP) { Z = MAXZINELP-1; }
      //G4cout << "Z= " << Z << G4endl;
      // Initialisation 
      if(!(data->GetElementData(Z))) { 
	Initialise(Z, dynParticle, path); 
      }
    }
    delete dynParticle;
  }
}

void G4ParticleInelasticXS::Initialise(G4int Z, G4DynamicParticle* dp, 
				       const char* p)
{
  if(data->GetElementData(Z)) { return; }
  const char* path = p;
  if(!p) {
    // check environment variable 
    // Build the complete string identifying the file with the data set
    path = getenv("G4PARTICLEXSDATA");
    if (!path) {
      G4Exception("G4ParticleInelasticXS::Initialise(..)","had013",
		  FatalException,
                  "Environment variable G4PARTICLEXSDATA is not defined");
      return;
    }
  }

  // upload element data 
  std::ostringstream ost;
  ost << path << "/" << particleName << "/inel" << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);
  /*
  G4cout << "G4ParticleInelasticXS::Initialise for Z= " << Z 
	 << " A= " << Amean << "  Amin= " << amin[Z] 
	 << "  Amax= " << amax[Z] << G4endl;
  */
  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << path << "/" << particleName << "/inel" << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v1); 
    }
  }

  // smooth transition 
  G4double sig1 = (*v)[v->GetVectorLength()-1];
  dp->SetKineticEnergy(v->GetMaxEnergy());
  G4double sig2 = 0.0;
  if(1 == Z) {
    if(fNucleon) {
      fNucleon->GetHadronNucleonXscNS(dp, proton);
      sig2 = fNucleon->GetInelasticHadronNucleonXsc();
    } else {
      nnXsection->GetZandACrossSection(dp, 1, 1);
      sig2 = nnXsection->GetInelasticGlauberGribovXsc();
    }
  } else {
    G4int Amean = G4lrint(fNist->GetAtomicMassAmu(Z));
    if(ggXsection) {
      ggXsection->GetIsoCrossSection(dp, Z, Amean);
      sig2 = ggXsection->GetInelasticGlauberGribovXsc();
    } else {
      nnXsection->GetZandACrossSection(dp, Z, Amean);
      sig2 = nnXsection->GetInelasticGlauberGribovXsc();
    }
  }
  if(sig2 > 0.) { coeff[Z] = sig1/sig2; } 
}

G4PhysicsVector* 
G4ParticleInelasticXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(warn) { 
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4ParticleInelasticXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4ParticleInelasticXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4ParticleInelasticXS::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  return v;
}
