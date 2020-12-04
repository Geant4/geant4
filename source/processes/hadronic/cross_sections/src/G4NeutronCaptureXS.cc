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
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4DynamicParticle.hh"
#include "G4ProductionCutsTable.hh"
#include "G4IsotopeList.hh"
#include "Randomize.hh"
#include "G4Log.hh" 

// factory
#include "G4CrossSectionFactory.hh"
//
G4_DECLARE_XS_FACTORY(G4NeutronCaptureXS);

G4ElementData* G4NeutronCaptureXS::data = nullptr;
G4String G4NeutronCaptureXS::gDataDirectory = "";

#ifdef G4MULTITHREADED
  G4Mutex G4NeutronCaptureXS::neutronCaptureXSMutex = G4MUTEX_INITIALIZER;
#endif

G4NeutronCaptureXS::G4NeutronCaptureXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   emax(20*CLHEP::MeV),elimit(1.0e-10*CLHEP::eV)
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronCaptureXS::G4NeutronCaptureXS: Initialise for Z < "
	    << MAXZCAPTURE << G4endl;
  }
  logElimit   = G4Log(elimit);
  isMaster    = false;
  temp.resize(13,0.0);
}

G4NeutronCaptureXS::~G4NeutronCaptureXS()
{
  if(isMaster) { delete data; data = nullptr; }
}

void G4NeutronCaptureXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronCaptureXS calculates the neutron capture cross sections\n"
          << "on nuclei using data from the high precision neutron database.\n"
          << "These data are simplified and smoothed over the resonance region\n"
          << "in order to reduce CPU time.  G4NeutronCaptureXS is set to zero\n"
          << "above 20 MeV for all targets. Cross section is zero also for Z>92.\n";
}
 
G4bool 
G4NeutronCaptureXS::IsElementApplicable(const G4DynamicParticle*, 
					G4int, const G4Material*)
{
  return true;
}

G4bool 
G4NeutronCaptureXS::IsIsoApplicable(const G4DynamicParticle*,
				    G4int, G4int,
				    const G4Element*, const G4Material*)
{
  return true;
}

G4double 
G4NeutronCaptureXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
					   G4int ZZ, const G4Material*)
{
  G4double xs = 0.0;
  G4double ekin = aParticle->GetKineticEnergy();
  if(ekin > emax) { return xs; }

  G4int Z = std::min(ZZ, MAXZCAPTURE-1);
  G4double logEkin = aParticle->GetLogKineticEnergy();
  if(ekin < elimit) { ekin = elimit; logEkin = logElimit; }

  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) { return xs; }

  const G4double e1 = pv->Energy(1);
  xs = (ekin >= e1) ? pv->LogVectorValue(ekin, logEkin) 
    : (*pv)[1]*std::sqrt(e1/ekin); 

#ifdef G4VERBOSE
  if(verboseLevel > 1){
    G4cout  << "Ekin= " << ekin/CLHEP::MeV 
            << " ElmXScap(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
  return xs;
}

G4double 
G4NeutronCaptureXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
				       G4int Z, G4int A,
				       const G4Isotope*, const G4Element*,
				       const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), 
                         aParticle->GetLogKineticEnergy(),
                         Z, A); 
}

G4double G4NeutronCaptureXS::IsoCrossSection(G4double eKin, G4double logE,
                                             G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  if(eKin > emax) { return xs; }

  G4int Z = std::min(ZZ, MAXZCAPTURE-1);
  G4double ekin = eKin;
  G4double logEkin = logE;
  if(ekin < elimit) { 
    ekin = elimit; 
    logEkin = logElimit; 
  }

  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) { return xs; }

  if(amin[Z] > 0 && A >= amin[Z] && A <= amax[Z]) {
    G4PhysicsVector* pviso = data->GetComponentDataByID(Z, A - amin[Z]);
    if(pviso != nullptr) { 
      const G4double e1 = pviso->Energy(1);
      xs = (ekin >= e1) ? pviso->LogVectorValue(ekin, logEkin)
	: (*pviso)[1]*std::sqrt(e1/ekin); 
#ifdef G4VERBOSE
      if(verboseLevel > 0) {
	G4cout << "G4NeutronCaptureXS::IsoXS: Ekin(MeV)= " << ekin/MeV 
	       << "  xs(b)= " << xs/barn 
	       << "  Z= " << Z << "  A= " << A << G4endl;
      }
#endif
      return xs;
    }
  }
  // isotope data are not available or applicable
  const G4double e1 = pv->Energy(1);
  xs = (ekin >= e1) ? pv->LogVectorValue(ekin, logEkin)
    : (*pv)[1]*std::sqrt(e1/ekin); 
#ifdef G4VERBOSE
  if(verboseLevel > 0) {
    G4cout << "G4NeutronCaptureXS::IsoXS: Ekin(MeV)= " << ekin/MeV 
           << "  xs(b)= " << xs/barn 
	   << "  Z= " << Z << "  A= " << A << " no iso XS" << G4endl;
  }
#endif
  return xs;
}

const G4Isotope* 
G4NeutronCaptureXS::SelectIsotope(const G4Element* anElement,
				  G4double kinEnergy, G4double logE)
{
  size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  // more than 1 isotope
  G4int Z = anElement->GetZasInt();

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;

  // is there isotope wise cross section?
  size_t j;
  if(0 == amin[Z] || Z >= MAXZCAPTURE) {
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

  if(nullptr == data) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&neutronCaptureXSMutex);
    if(nullptr == data) { 
#endif
      isMaster = true;
      data = new G4ElementData(); 
      data->SetName("NeutronCapture");
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&neutronCaptureXSMutex);
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
	G4int Z = std::max(1,std::min(((*elmVec)[ie])->GetZasInt(), MAXZCAPTURE-1));
	if(nullptr == data->GetElementData(Z)) { Initialise(Z); }
      }
    }
  }
}

const G4String& G4NeutronCaptureXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    char* path = std::getenv("G4PARTICLEXSDATA");
    if (path) {
      std::ostringstream ost;
      ost << path << "/neutron/cap";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4NeutronCaptureXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEXSDATA is not defined");
    }
  }
  return gDataDirectory;
}

void G4NeutronCaptureXS::InitialiseOnFly(G4int Z)
{
#ifdef G4MULTITHREADED
   G4MUTEXLOCK(&neutronCaptureXSMutex);
   if(nullptr == data->GetElementData(Z)) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&neutronCaptureXSMutex);
#endif
}

void G4NeutronCaptureXS::Initialise(G4int Z)
{
  if(nullptr != data->GetElementData(Z)) { return; }

  // upload element data 
  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);

  // upload isotope data
  if(amin[Z] > 0) {
    size_t nmax = (size_t)(amax[Z]-amin[Z]+1);
    data->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << gDataDirectory << Z << "_" << A;
      v = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v); 
    }
  }
}
 
G4PhysicsVector* 
G4NeutronCaptureXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    if(warn) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4NeutronCaptureXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
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
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  return v;
}
