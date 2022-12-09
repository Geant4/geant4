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
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"
#include "G4AutoLock.hh"

#include <fstream>
#include <sstream>

G4PhysicsVector* G4NeutronElasticXS::data[] = {nullptr};
G4double G4NeutronElasticXS::coeff[] = {0.0};
G4String G4NeutronElasticXS::gDataDirectory = "";

namespace
{
  G4Mutex nElasticXSMutex = G4MUTEX_INITIALIZER;
}

G4NeutronElasticXS::G4NeutronElasticXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   neutron(G4Neutron::Neutron())
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4NeutronElasticXS::G4NeutronElasticXS Initialise for Z < " 
	    << MAXZEL << G4endl;
  }
  ggXsection = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if(ggXsection == nullptr) ggXsection = new G4ComponentGGHadronNucleusXsc();
  SetForAllAtomsAndEnergies(true);
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
  return false;
}

G4double 
G4NeutronElasticXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
					   G4int Z, const G4Material*)
{
  return ElementCrossSection(aParticle->GetKineticEnergy(), aParticle->GetLogKineticEnergy(), Z);
}

G4double
G4NeutronElasticXS::ComputeCrossSectionPerElement(G4double ekin, G4double loge,
						  const G4ParticleDefinition*,
						  const G4Element* elm,
						  const G4Material*)
{
  return ElementCrossSection(ekin, loge, elm->GetZasInt());
}

G4double G4NeutronElasticXS::ElementCrossSection(G4double ekin, G4double loge, G4int ZZ)
{
  G4int Z = (ZZ >= MAXZEL) ? MAXZEL - 1 : ZZ;
  auto pv = GetPhysicsVector(Z);

  G4double xs = (ekin <= pv->GetMaxEnergy()) ? pv->LogVectorValue(ekin, loge) 
    : coeff[Z]*ggXsection->GetElasticElementCrossSection(neutron, ekin,
                                                         Z, aeff[Z]);

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ",  nElmXSel(b)= " << xs/CLHEP::barn 
	    << G4endl;
  }
#endif
  return xs;
}

G4double
G4NeutronElasticXS::ComputeIsoCrossSection(G4double ekin, G4double loge,
				           const G4ParticleDefinition*,
				           G4int Z, G4int A,
				           const G4Isotope*, const G4Element*,
				           const G4Material*)
{
  return ElementCrossSection(ekin, loge, Z)*A/aeff[Z];
}

G4double
G4NeutronElasticXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
				       G4int Z, G4int A,
				       const G4Isotope*, const G4Element*,
				       const G4Material*)
{
  return ElementCrossSection(aParticle->GetKineticEnergy(),
			     aParticle->GetLogKineticEnergy(), Z)*A/aeff[Z];

}

const G4Isotope* G4NeutronElasticXS::SelectIsotope(
      const G4Element* anElement, G4double, G4double)
{
  G4int nIso = (G4int)anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;

  // isotope wise cross section not used
  for (G4int j=0; j<nIso; ++j) {
    sum += abundVector[j];
    if(q <= sum) {
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
    G4AutoLock l(&nElasticXSMutex);
    if(0. == coeff[0]) { 
      coeff[0] = 1.0;
      isMaster = true;
      FindDirectoryPath();
    }
    l.unlock();
  }

  // it is possible re-initialisation for the second run
  if(isMaster) {

    // Access to elements
    const G4ElementTable* table = G4Element::GetElementTable();
    for ( auto & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZEL-1) );
      if ( nullptr == data[Z] ) { Initialise(Z); }
    }
  }
}

const G4String& G4NeutronElasticXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4PARTICLEXSDATA");
    if (nullptr != path) {
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
  if(nullptr == data[Z]) {
    G4AutoLock l(&nElasticXSMutex);
    if(nullptr == data[Z]) { Initialise(Z); }
    l.unlock();
  }
}

void G4NeutronElasticXS::Initialise(G4int Z)
{
  if(data[Z] != nullptr) { return; }

  // upload data from file
  data[Z] = new G4PhysicsLogVector();

  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  std::ifstream filein(ost.str().c_str());
  if (!filein.is_open()) {
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
  G4double sig2  = ggXsection->GetElasticElementCrossSection(neutron, 
                               ehigh, Z, aeff[Z]);
  coeff[Z] = (sig2 > 0.) ? sig1/sig2 : 1.0;  
}
