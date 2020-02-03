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
// File name:    G4NeutronElasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
//

#include "G4NeutronElasticXS.hh"
#include "G4Neutron.hh"
#include "G4DynamicParticle.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsVector.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <sstream>

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4NeutronElasticXS);

using namespace std;

G4PhysicsVector* G4NeutronElasticXS::data[] = {nullptr};
G4double G4NeutronElasticXS::coeff[] = {0.0};
G4double G4NeutronElasticXS::aeff[]  = {1.0};
G4String G4NeutronElasticXS::gDataDirectory = "";

#ifdef G4MULTITHREADED
  G4Mutex G4NeutronElasticXS::neutronElasticXSMutex = G4MUTEX_INITIALIZER;
#endif

G4NeutronElasticXS::G4NeutronElasticXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   ggXsection(nullptr),
   neutron(G4Neutron::Neutron()),
   isMaster(false)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronElasticXS::G4NeutronElasticXS Initialise for Z < " 
	    << MAXZEL << G4endl;
  }
  nist = G4NistManager::Instance();
  ggXsection = new G4ComponentGGHadronNucleusXsc();
  SetForAllAtomsAndEnergies(true);
  temp.resize(13,0.0);
}

G4NeutronElasticXS::~G4NeutronElasticXS()
{
  if(isMaster) {
    for(G4int i=0; i<MAXZEL; ++i) {
      delete data[i];
      data[i] = nullptr;
    }
  }
}

void G4NeutronElasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronElasticXS calculates the neutron elastic scattering\n"
          << "cross section on nuclei using data from the high precision\n"
          << "neutron database.  These data are simplified and smoothed over\n"
          << "the resonance region in order to reduce CPU time.\n"
          << "For high energies Glauber-Gribiv cross section is used.\n";
}

G4bool 
G4NeutronElasticXS::IsElementApplicable(const G4DynamicParticle*, 
					G4int, const G4Material*)
{
  return true;
}

G4bool G4NeutronElasticXS::IsIsoApplicable(const G4DynamicParticle*,
                                           G4int, G4int,
                                           const G4Element*, const G4Material*)
{
  return true;
}

G4double 
G4NeutronElasticXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
					   G4int ZZ, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();

  G4int Z = (ZZ >= MAXZEL) ? MAXZEL - 1 : ZZ; 

  auto pv = GetPhysicsVector(Z);
  if(!pv) { return xs; }
  //  G4cout  << "G4NeutronElasticXS::GetCrossSection e= " << ekin 
  // << " Z= " << Z << G4endl;

  if(ekin <= pv->Energy(0)) { 
    xs = (*pv)[0];
  } else if(ekin <= pv->GetMaxEnergy()) { 
    xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy()); 
  } else {          
    xs = coeff[Z]*ggXsection->GetElasticElementCrossSection(neutron, 
                  ekin, Z, aeff[Z]);
  }

  if(verboseLevel > 1) {
    G4cout  << "Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ",  nElmXSel(b)= " << xs/CLHEP::barn 
	    << G4endl;
  }
  return xs;
}

G4double G4NeutronElasticXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int Z, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), 
                         aParticle->GetLogKineticEnergy(), Z, A);
}

G4double 
G4NeutronElasticXS::IsoCrossSection(G4double ekin, G4double logekin, 
                                    G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  G4int Z = (ZZ >= MAXZEL) ? MAXZEL - 1 : ZZ; 

  // tritium and He3 
  if(3 == A) {
    return ggXsection->GetElasticElementCrossSection(neutron, ekin, Z, A);
  }
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  auto pv = GetPhysicsVector(Z);
  if(!pv) { return xs; }

  if(ekin <= pv->Energy(0)) { 
    xs = (*pv)[0];
  } else if(ekin <= pv->GetMaxEnergy()) { 
    xs = pv->LogVectorValue(ekin, logekin); 
  } else {          
    xs = coeff[Z]*ggXsection->GetElasticElementCrossSection(neutron, 
                  ekin, Z, aeff[Z]);
  }
  xs *= A/aeff[Z];
  if(verboseLevel > 1) {
    G4cout  << "G4NeutronElasticXS::IsoXS: Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
  return xs;
}

const G4Isotope* G4NeutronElasticXS::SelectIsotope(
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

  // isotope wise cross section not used
  if(anElement->GetNaturalAbundanceFlag()) {
    for (j=0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = anElement->GetIsotope(j);
	break;
      }
    }
    return iso;
  }

  // use isotope cross sections
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
  for (j = 0; j<nIso; ++j) {
    if(temp[j] >= sum) {
      iso = anElement->GetIsotope(j);
      break;
    }
  }
  return iso;
}

void 
G4NeutronElasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0){
    G4cout << "G4NeutronElasticXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "neutron") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronElasticXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }
  if(0. == coeff[0]) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&neutronElasticXSMutex);
    if(0. == coeff[0]) { 
#endif
      coeff[0] = 1.0;
      isMaster = true;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&neutronElasticXSMutex);
#endif
  }

  // it is possible re-initialisation for the second run
  if(isMaster) {

    // Access to elements
    auto theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();
    for(size_t j=0; j<numOfCouples; ++j) {
      auto mat = theCoupleTable->GetMaterialCutsCouple(j)->GetMaterial();
      auto elmVec = mat->GetElementVector();
      size_t numOfElem = mat->GetNumberOfElements();
      for (size_t ie = 0; ie < numOfElem; ++ie) {
	G4int Z = std::max(1,std::min(((*elmVec)[ie])->GetZasInt(), MAXZEL-1));
	if(!data[Z]) { Initialise(Z); }
      }
    }
  }
}

G4PhysicsVector* G4NeutronElasticXS::GetPhysicsVector(G4int Z)
{
  if(!data[Z]) { InitialiseOnFly(Z); }
  return data[Z];
}

const G4String& G4NeutronElasticXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    char* path = std::getenv("G4PARTICLEXSDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/neutron/el";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4NeutronElasticXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEXSDATA is not defined");
    }
  }
  return gDataDirectory;
}

void G4NeutronElasticXS::InitialiseOnFly(G4int Z)
{
#ifdef G4MULTITHREADED
   G4MUTEXLOCK(&neutronElasticXSMutex);
   if(!data[Z]) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&neutronElasticXSMutex);
#endif
}

void G4NeutronElasticXS::Initialise(G4int Z)
{
  if(data[Z]) { return; }

  // upload data from file
  data[Z] = new G4PhysicsLogVector();

  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not opened!";
    G4Exception("G4NeutronElasticXS::Initialise(..)","had014",
                FatalException, ed, "Check G4PARTICLEXSDATA");
    return;
  }
  if(verboseLevel > 1) {
    G4cout << "file " << ost.str() 
	   << " is opened by G4NeutronElasticXS" << G4endl;
  }
    
  // retrieve data from DB
  if(!data[Z]->Retrieve(filein, true)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not retrieved!";
    G4Exception("G4NeutronElasticXS::Initialise(..)","had015",
		FatalException, ed, "Check G4PARTICLEXSDATA");
    return;
  }
  // smooth transition 
  G4double sig1  = (*(data[Z]))[data[Z]->GetVectorLength()-1];
  G4double ehigh = data[Z]->GetMaxEnergy();
  aeff[Z] = nist->GetAtomicMassAmu(Z);
  G4double sig2  = ggXsection->GetElasticElementCrossSection(neutron, 
                               ehigh, Z, aeff[Z]);
  if(sig2 > 0.) { coeff[Z] = sig1/sig2; } 
}
