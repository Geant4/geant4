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
// File name:    G4GammaNuclearXS
//
// Author  V.Ivantchenko, Geant4, 22 October 2020
//
// Modifications:
//

#include "G4GammaNuclearXS.hh"
#include "G4Gamma.hh"
#include "G4DynamicParticle.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"

#include <fstream>
#include <sstream>

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4GammaNuclearXS);

G4PhysicsVector* G4GammaNuclearXS::data[] = {nullptr};
G4double G4GammaNuclearXS::coeff[] = {0.0};
G4String G4GammaNuclearXS::gDataDirectory = "";

#ifdef G4MULTITHREADED
  G4Mutex G4GammaNuclearXS::gNuclearXSMutex = G4MUTEX_INITIALIZER;
#endif

G4GammaNuclearXS::G4GammaNuclearXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   ggXsection(nullptr),
   gamma(G4Gamma::Gamma()),
   isMaster(false)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4GammaNuclearXS::G4GammaNuclearXS Initialise for Z < " 
	    << MAXZGAMMAN << G4endl;
  }
  ggXsection = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet("PhotoNuclearXS");
  if(ggXsection == nullptr) ggXsection = new G4PhotoNuclearCrossSection();
  SetForAllAtomsAndEnergies(true);
}

G4GammaNuclearXS::~G4GammaNuclearXS()
{
  if(isMaster) {
    for(G4int i=0; i<MAXZGAMMAN; ++i) {
      delete data[i];
      data[i] = nullptr;
    }
  }
}

void G4GammaNuclearXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4GammaNuclearXS calculates the gamma nuclear\n"
          << "cross section on nuclei using data from the high precision\n"
          << "LEND gamma database. The data are simplified and smoothed over\n"
          << "the resonance region in order to reduce CPU time.\n"
          << "For high energies Glauber-Gribiv cross section is used.\n";
}

G4bool 
G4GammaNuclearXS::IsElementApplicable(const G4DynamicParticle*, 
                                      G4int, const G4Material*)
{
  return true;
}

G4bool G4GammaNuclearXS::IsIsoApplicable(const G4DynamicParticle*,
                                         G4int, G4int,
                                         const G4Element*, const G4Material*)
{
  return true;
}

G4double 
G4GammaNuclearXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
                                         G4int ZZ, const G4Material* mat)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();

  G4int Z = (ZZ >= MAXZGAMMAN) ? MAXZGAMMAN - 1 : ZZ; 

  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) { return xs; }
  //  G4cout  << "G4GammaNuclearXS::GetCrossSection e= " << ekin 
  // << " Z= " << Z << G4endl;

  if(ekin <= pv->GetMaxEnergy()) { 
    xs = pv->LogVectorValue(ekin, aParticle->GetLogKineticEnergy()); 
  } else {          
    xs = coeff[Z]*ggXsection->GetElementCrossSection(aParticle, Z, mat);
  }

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ",  nElmXS(b)= " << xs/CLHEP::barn 
	    << G4endl;
  }
#endif
  return xs;
}

G4double G4GammaNuclearXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle, 
	 G4int Z, G4int A,
	 const G4Isotope*, const G4Element*,
	 const G4Material* mat)
{
  return GetElementCrossSection(aParticle, Z, mat) * A/aeff[Z];
}

const G4Isotope* G4GammaNuclearXS::SelectIsotope(
      const G4Element* anElement, G4double, G4double)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;
  size_t j;

  // isotope wise cross section not used
  for (j=0; j<nIso; ++j) {
    sum += abundVector[j];
    if(q <= sum) {
      iso = anElement->GetIsotope(j);
      break;
    }
  }
  return iso;
}

void 
G4GammaNuclearXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0){
    G4cout << "G4GammaNuclearXS::BuildPhysicsTable for " 
	   << p.GetParticleName() << G4endl;
  }
  if(p.GetParticleName() != "gamma") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only gamma is allowed";
    G4Exception("G4GammaNuclearXS::BuildPhysicsTable(..)","had012",
		FatalException, ed, "");
    return; 
  }
  if(0. == coeff[0]) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&gNuclearXSMutex);
    if(0. == coeff[0]) { 
#endif
      coeff[0] = 1.0;
      isMaster = true;
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&gNuclearXSMutex);
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
	G4int Z = std::max(1,std::min(((*elmVec)[ie])->GetZasInt(), MAXZGAMMAN-1));
	if(data[Z] == nullptr) { Initialise(Z); }
      }
    }
  }
}

const G4String& G4GammaNuclearXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    char* path = std::getenv("G4PARTICLEXSDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/gamma/inel";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4GammaNuclearXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEXSDATA is not defined");
    }
  }
  return gDataDirectory;
}

void G4GammaNuclearXS::InitialiseOnFly(G4int Z)
{
#ifdef G4MULTITHREADED
   G4MUTEXLOCK(&gNuclearXSMutex);
   if(data[Z] == nullptr) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&gNuclearXSMutex);
#endif
}

void G4GammaNuclearXS::Initialise(G4int Z)
{
  if(data[Z] != nullptr) { return; }

  // upload data from file
  data[Z] = new G4PhysicsLogVector();

  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not opened!";
    G4Exception("G4GammaNuclearXS::Initialise(..)","had014",
                FatalException, ed, "Check G4PARTICLEXSDATA");
    return;
  }
  if(verboseLevel > 1) {
    G4cout << "file " << ost.str() 
	   << " is opened by G4GammaNuclearXS" << G4endl;
  }
    
  // retrieve data from DB
  if(!data[Z]->Retrieve(filein, true)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not retrieved!";
    G4Exception("G4GammaNuclearXS::Initialise(..)","had015",
		FatalException, ed, "Check G4PARTICLEXSDATA");
    return;
  }
  // smooth transition 
  G4Material* mat(nullptr);
  G4ThreeVector mom(0.0,0.0,1.0);
  G4DynamicParticle dp(gamma, mom, data[Z]->GetMaxEnergy());
  G4double sig1  = (*(data[Z]))[data[Z]->GetVectorLength()-1];
  G4double sig2  = ggXsection->GetElementCrossSection(&dp, Z, mat);
  coeff[Z] = (sig2 > 0.) ? sig1/sig2 : 1.0; 
}
