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
#include "G4HadronicParameters.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"
#include "G4AutoLock.hh"

#include <fstream>
#include <sstream>

G4ElementData* G4NeutronElasticXS::data = nullptr;
G4ElementData* G4NeutronElasticXS::dataR = nullptr;
G4double G4NeutronElasticXS::coeff[] = {1.0};
G4String G4NeutronElasticXS::gDataDirectory = "";

static std::once_flag applyOnce;

namespace
{
  G4Mutex nElasticXSMutex = G4MUTEX_INITIALIZER;
}

G4NeutronElasticXS::G4NeutronElasticXS() 
 : G4VCrossSectionDataSet(Default_Name()),
   neutron(G4Neutron::Neutron())
{
  if (verboseLevel > 0) {
    G4cout  << "G4NeutronElasticXS::G4NeutronElasticXS Initialise for Z < " 
	    << MAXZEL << G4endl;
  }
  ggXsection = 
    G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if (ggXsection == nullptr)
    ggXsection = new G4ComponentGGHadronNucleusXsc();
  SetForAllAtomsAndEnergies(true);
  if (nullptr == data) { 
    data = new G4ElementData(MAXZEL);
    data->SetName("nElastic");
    dataR = new G4ElementData(MAXZEL);
    dataR->SetName("nRElastic");
    FindDirectoryPath();
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
  return ElementCrossSection(aParticle->GetKineticEnergy(),
			     aParticle->GetLogKineticEnergy(), Z);
}

G4double
G4NeutronElasticXS::ComputeCrossSectionPerElement(G4double ekin, G4double loge,
						  const G4ParticleDefinition*,
						  const G4Element* elm,
						  const G4Material*)
{
  return ElementCrossSection(ekin, loge, elm->GetZasInt());
}

G4double
G4NeutronElasticXS::ElementCrossSection(G4double ekin, G4double loge, G4int ZZ)
{
  G4int Z = std::min(ZZ, MAXZEL-1);
  G4double xs;
  G4bool done{false};
  
  // data from the resonance region
  if (fRfilesEnabled) {
    auto pv = GetPhysicsVectorR(Z);
    if (nullptr != pv && ekin < el_max_r_e[Z]) {
      xs = pv->LogVectorValue(ekin, loge);
      done = true;
    }
  }
  // data above the resonance region
  if (!done) { 
    auto pv = GetPhysicsVector(Z);
    xs = (ekin <= pv->GetMaxEnergy()) ? pv->LogVectorValue(ekin, loge)
      : coeff[Z]*ggXsection->GetElasticElementCrossSection(neutron, ekin,
							   Z, aeff[Z]);
  }

#ifdef G4VERBOSE
  if (verboseLevel > 1) {
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
				           G4int ZZ, G4int A,
				           const G4Isotope*, const G4Element*,
				           const G4Material*)
{
  G4int Z = std::min(ZZ, MAXZEL-1);
  return ElementCrossSection(ekin, loge, Z)*A/aeff[Z];
}

G4double
G4NeutronElasticXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
				       G4int ZZ, G4int A,
				       const G4Isotope*, const G4Element*,
				       const G4Material*)
{
  G4int Z = std::min(ZZ, MAXZEL-1);
  return ElementCrossSection(aParticle->GetKineticEnergy(),
			     aParticle->GetLogKineticEnergy(), Z)*A/aeff[Z];

}

const G4Isotope* G4NeutronElasticXS::SelectIsotope(
      const G4Element* anElement, G4double, G4double)
{
  G4int nIso = (G4int)anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

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

  fRfilesEnabled = G4HadronicParameters::Instance()->UseRFilesForXS();

  // initialise static tables only once
  std::call_once(applyOnce, [this]() { isInitializer = true; });

  if (isInitializer) {
    G4AutoLock l(&nElasticXSMutex);

    // Access to elements
    const G4ElementTable* table = G4Element::GetElementTable();
    for ( auto const & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZEL-1) );
      if ( nullptr == data->GetElementData(Z) ) { Initialise(Z); }
    }
    l.unlock();
  }
}

const G4String& G4NeutronElasticXS::FindDirectoryPath()
{
  // build the complete string identifying the file with the data set
  if (gDataDirectory.empty()) {
    std::ostringstream ost;
    ost << G4HadronicParameters::Instance()->GetDirPARTICLEXS() << "/neutron/";
    gDataDirectory = ost.str();
  }
  return gDataDirectory;
}

void G4NeutronElasticXS::InitialiseOnFly(G4int Z)
{
  G4AutoLock l(&nElasticXSMutex);
  Initialise(Z);
  l.unlock();
}

void G4NeutronElasticXS::Initialise(G4int Z)
{
  if (nullptr != data->GetElementData(Z)) { return; }

  // upload element data
  std::ostringstream ost;
  ost << FindDirectoryPath() << "el" << Z;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);

  G4PhysicsVector* vr = nullptr;
  if (fRfilesEnabled) {
    std::ostringstream ostr;
    ostr << FindDirectoryPath() << "Rel" << Z;
    vr = RetrieveVector(ostr, false);
    dataR->InitialiseForElement(Z, vr);
  }
  
  // smooth transition 
  G4double sig1 = (*v)[v->GetVectorLength()-1];
  G4double ehigh = v->GetMaxEnergy();
  G4double sig2 =
    ggXsection->GetElasticElementCrossSection(neutron, ehigh, Z, aeff[Z]);
  coeff[Z] = (sig2 > 0.) ? sig1/sig2 : 1.0;  
}

G4PhysicsVector* 
G4NeutronElasticXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!filein.is_open()) {
    if (warn) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not opened!";
      G4Exception("G4NeutronElasticXS::RetrieveVector(..)","had014",
                  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  } else {
    if (verboseLevel > 1) {
      G4cout << "File " << ost.str() 
             << " is opened by G4NeutronElasticXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if (!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not retrieved!";
      G4Exception("G4NeutronElasticXS::RetrieveVector(..)","had015",
                  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  return v;
}
