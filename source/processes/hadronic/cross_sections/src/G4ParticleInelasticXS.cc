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
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4Proton.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"
#include "G4HadronicParameters.hh"

#include <fstream>
#include <sstream>

G4ElementData* G4ParticleInelasticXS::data[] = {nullptr};
G4double G4ParticleInelasticXS::coeff[MAXZINELP][5] = {{1.0}, {1.0}, {1.0}, {1.0}, {1.0}};
G4String G4ParticleInelasticXS::gDataDirectory[] = {""};

#ifdef G4MULTITHREADED
  G4Mutex G4ParticleInelasticXS::particleInelasticXSMutex = G4MUTEX_INITIALIZER;
#endif

G4ParticleInelasticXS::G4ParticleInelasticXS(const G4ParticleDefinition* part) 
  : G4VCrossSectionDataSet("G4ParticleInelasticXS"),
    highEnergyXsection(nullptr),
    particle(part),
    index(0),
    isMaster(false)
{
  if(!part) {
    G4Exception("G4ParticleInelasticXS::G4ParticleInelasticXS(..)","had015",
		FatalException, "NO particle definition in constructor");
  } else {
    verboseLevel = 0;
    const G4String& particleName = particle->GetParticleName();
    if(verboseLevel > 1) {
      G4cout << "G4ParticleInelasticXS::G4ParticleInelasticXS for " 
	     << particleName << " on atoms with Z < " << MAXZINELP << G4endl;
    }
    if(particleName == "proton") {
      highEnergyXsection = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
      if(highEnergyXsection == nullptr) {
	highEnergyXsection = new G4ComponentGGHadronNucleusXsc();
      }
    } else {
      highEnergyXsection = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov Nucl-nucl");
      if(highEnergyXsection == nullptr) {
	highEnergyXsection = new G4ComponentGGNuclNuclXsc();
      }
      if(particleName == "deuteron") index = 1; 
      else if(particleName == "triton") index = 2; 
      else if(particleName == "He3") index = 3; 
      else if(particleName == "alpha") index = 4;
      else {
        G4ExceptionDescription ed;
        ed << particleName << " is a wrong particle type";
	G4Exception("G4ParticleInelasticXS::BuildPhysicsTable(..)","had012",
		    FatalException, ed, "");
      } 
    }
  }
  SetForAllAtomsAndEnergies(true);
  temp.resize(13,0.0);
}

G4ParticleInelasticXS::~G4ParticleInelasticXS()
{
  if(isMaster) { 
    delete data[index]; 
    data[index] = nullptr;
  }
}

void G4ParticleInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4ParticleInelasticXS calculates n, p, d, t, he3, he4 inelastic\n"
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

  auto pv = GetPhysicsVector(Z);
  if(nullptr == pv) { return xs; }
  //  G4cout  << "G4ParticleInelasticXS::GetCrossSection e= " << ekin 
  //  << " Z= " << Z << G4endl;

  if(ekin <= pv->GetMaxEnergy()) { 
    xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy()); 
  } else {
    xs = coeff[Z][index]*highEnergyXsection->GetInelasticElementCrossSection(particle,
			        ekin, Z, aeff[Z]);
  }

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "ElmXS: Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << " xs(bn)= " << xs/CLHEP::barn << " element data for "
	    << particle->GetParticleName()
            << " idx= " << index << G4endl;
  }
#endif
  return xs;
}

G4double G4ParticleInelasticXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int Z, G4int A,
	 const G4Isotope*, const G4Element*, const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), 
                         aParticle->GetLogKineticEnergy(),Z, A);
}

G4double 
G4ParticleInelasticXS::IsoCrossSection(G4double ekin, G4double logE,
                                       G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  G4int Z = (ZZ >= MAXZINELP) ? MAXZINELP - 1 : ZZ; 
  /*  
  G4cout << "G4ParticleInelasticXS: IsoCrossSection  Z= " 
         << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) { return xs; }

  // compute isotope cross section if applicable 
  const G4double emax = pv->GetMaxEnergy(); 
  if(ekin <= emax && amin[Z]>0 && A >= amin[Z] && A <= amax[Z]) {
    auto pviso = data[index]->GetComponentDataByIndex(Z, A - amin[Z]);
    if(pviso != nullptr) { 
      xs = pviso->LogVectorValue(ekin, logE);
#ifdef G4VERBOSE
      if(verboseLevel > 1) {
	G4cout << "G4ParticleInelasticXS::IsoXS: for " 
               << particle->GetParticleName() << " Ekin(MeV)= " 
               << ekin/CLHEP::MeV << "  xs(b)= " << xs/CLHEP::barn 
	       << "  Z= " << Z << "  A= " << A 
               << " idx= " << index << G4endl;
      }
#endif
      return xs;
    }
  }
  // use element x-section
  if(ekin <= emax) { 
    xs = pv->LogVectorValue(ekin, logE); 
  } else {
    xs = coeff[Z][index] *
      highEnergyXsection->GetInelasticElementCrossSection(particle,
			  ekin, Z, aeff[Z]);
  }
  xs *= A/aeff[Z];
#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "IsoXS for " << particle->GetParticleName() 
	    << " Target Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << " xs(bn)= " << xs/CLHEP::barn
            << " idx= " << index << G4endl;
  }
#endif
  return xs;
}

const G4Isotope* G4ParticleInelasticXS::SelectIsotope(
		 const G4Element* anElement, G4double kinEnergy, G4double logE)
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
  size_t j;

  // isotope wise cross section not available
  if(0 == amin[Z] || Z >= MAXZINELP) {
    for (j=0; j<nIso; ++j) {
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
    sum += abundVector[j]*IsoCrossSection(kinEnergy, logE, Z, 
					  anElement->GetIsotope(j)->GetN());
    temp[j] = sum;
  }
  sum *= q;
  for (j=0; j<nIso; ++j) {
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
       << particle->GetParticleName() << " is expected";
    G4Exception("G4ParticleInelasticXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }

  G4int fact = (p.GetParticleName() == "proton") ? 1 : 256;
  SetMaxKinEnergy(G4HadronicParameters::Instance()->GetMaxEnergy() * fact);

  if(data[index] == nullptr) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&particleInelasticXSMutex);
    if(data[index] == nullptr) { 
#endif
      isMaster = true;
      data[index] = new G4ElementData(); 
      data[index]->SetName(particle->GetParticleName() + "Inelastic");
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&particleInelasticXSMutex);
#endif
  }

  // it is possible re-initialisation for the new run
  if(isMaster) {

    // Access to elements
    auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();
    for(size_t j=0; j<numOfCouples; ++j) {
      auto mat = theCoupleTable->GetMaterialCutsCouple(j)->GetMaterial();
      auto elmVec = mat->GetElementVector();
      size_t numOfElem = mat->GetNumberOfElements();
      for (size_t ie = 0; ie < numOfElem; ++ie) {
	G4int Z = std::max(1,std::min(((*elmVec)[ie])->GetZasInt(), MAXZINELP-1));
	if(nullptr == data[index]->GetElementData(Z)) { Initialise(Z); }
      }
    }
  }
}

const G4String& G4ParticleInelasticXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory[index].empty()) {
    char* path = std::getenv("G4PARTICLEXSDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/" << particle->GetParticleName() << "/inel";
      gDataDirectory[index] = ost.str();
    } else {
      G4Exception("G4NeutronInelasticXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEXSDATA is not defined");
    }
  }
  return gDataDirectory[index];
}

void G4ParticleInelasticXS::InitialiseOnFly(G4int Z)
{
#ifdef G4MULTITHREADED
   G4MUTEXLOCK(&particleInelasticXSMutex);
   if(nullptr == data[index]->GetElementData(Z)) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&particleInelasticXSMutex);
#endif
}

void G4ParticleInelasticXS::Initialise(G4int Z)
{
  if(nullptr != data[index]->GetElementData(Z)) { return; }

  // upload element data 
  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data[index]->InitialiseForElement(Z, v);
  /*
  G4cout << "G4ParticleInelasticXS::Initialise for Z= " << Z 
	 << " idx= " << index
	 << "  Amin= " << amin[Z] 
	 << "  Amax= " << amax[Z] 
	 << "  " << FindDirectoryPath() << G4endl;
  */
  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data[index]->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << FindDirectoryPath() << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      data[index]->AddComponent(Z, A, v1); 
      /*
      G4cout << "  Isotope x-section Z= " << Z << " A= " << A 
	     << " v1= " << v1 << G4endl;
      */
    }
  }
  // smooth transition 
  G4double sig1  = (*v)[v->GetVectorLength()-1];
  G4double ehigh = v->GetMaxEnergy();
  G4double sig2 = highEnergyXsection->GetInelasticElementCrossSection(
                  particle, ehigh, Z, aeff[Z]);
  coeff[Z][index] = (sig2 > 0.) ? sig1/sig2 : 1.0;
  /*
  G4cout << "G4ParticleInelasticXS: index= " << index 
         << " Z= " << Z << " coeff= " << coeff[Z][index] << G4endl; 
  */
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
