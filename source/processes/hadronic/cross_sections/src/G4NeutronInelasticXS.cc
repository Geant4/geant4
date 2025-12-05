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
// File name:    G4NeutronInelasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//

#include "G4NeutronInelasticXS.hh"
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
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"
#include "G4NuclearRadii.hh"
#include "G4AutoLock.hh"

#include <fstream>
#include <sstream>
#include <thread>

G4double G4NeutronInelasticXS::coeff[] = {1.0};
G4double G4NeutronInelasticXS::lowcoeff[] = {1.0};
G4ElementData* G4NeutronInelasticXS::data = nullptr;
G4ElementData* G4NeutronInelasticXS::dataR = nullptr;
G4String G4NeutronInelasticXS::gDataDirectory = "";

static std::once_flag applyOnce;

namespace
{
  G4Mutex nInelasticXSMutex = G4MUTEX_INITIALIZER;
}

G4NeutronInelasticXS::G4NeutronInelasticXS() 
  : G4VCrossSectionDataSet(Default_Name()),
    neutron(G4Neutron::Neutron()),
    lowElimit(1.0e-7*CLHEP::eV)
{
  verboseLevel = 0;
  if (verboseLevel > 0) {
    G4cout << "G4NeutronInelasticXS::G4NeutronInelasticXS Initialise for Z < " 
            << MAXZINEL << G4endl;
  }
  loglowElimit = G4Log(lowElimit);
  ggXsection =
    G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if (ggXsection == nullptr)
    ggXsection = new G4ComponentGGHadronNucleusXsc();

  if (nullptr == data) { 
    data = new G4ElementData(MAXZINEL);
    data->SetName("nInelastic");
    dataR = new G4ElementData(MAXZINEL);
    dataR->SetName("nRInelastic");
    FindDirectoryPath();
    for (G4int Z=1; Z<MAXZINEL; ++Z) { Initialise(Z); }
  }

  SetForAllAtomsAndEnergies(true);
}

void G4NeutronInelasticXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4NeutronInelasticXS calculates the neutron inelastic scattering\n"
          << "cross section on nuclei using data from the high precision\n"
          << "neutron database.  These data are simplified and smoothed over\n"
          << "the resonance region in order to reduce CPU time.\n"
          << "For high energy Glauber-Gribov cross section model is used\n";
}

G4bool 
G4NeutronInelasticXS::IsElementApplicable(const G4DynamicParticle*, 
                                          G4int, const G4Material*)
{
  return true;
}

G4bool 
G4NeutronInelasticXS::IsIsoApplicable(const G4DynamicParticle*,
                                      G4int, G4int,
                                      const G4Element*, const G4Material*)
{
  return true;
}

G4double 
G4NeutronInelasticXS::GetElementCrossSection(const G4DynamicParticle* aParticle,
                                             G4int Z, const G4Material*)
{
  return ElementCrossSection(aParticle->GetKineticEnergy(),
                             aParticle->GetLogKineticEnergy(), Z);
}

G4double
G4NeutronInelasticXS::ComputeCrossSectionPerElement(G4double ekin, G4double loge,
                                                    const G4ParticleDefinition*,
                                                    const G4Element* elm,
                                                    const G4Material*)
{
  return ElementCrossSection(ekin, loge, elm->GetZasInt());
}

G4double
G4NeutronInelasticXS::ElementCrossSection(G4double eKin, G4double logE, G4int ZZ)
{
  G4int Z = std::min(ZZ, MAXZINEL-1);
  G4double ekin = eKin;
  G4double loge = logE;
  G4double xs{0.0};
  G4bool done{false};

  // very low energy limit
  if (ekin < lowElimit) { 
    ekin = lowElimit; 
    loge = loglowElimit; 
  }

  // data from the resonance region
  if (fRfilesEnabled) {
    auto pv = GetPhysicsVectorR(Z);
    if (nullptr != pv && ekin < inel_max_r_e[Z]) {
      xs = pv->LogVectorValue(ekin, loge);
      done = true;
    }
  }
  // data above the resonance region pv should always defined
  if (!done) {
    auto pv = GetPhysicsVector(Z);
    if (ekin <= pv->GetMaxEnergy()) {
      xs = pv->LogVectorValue(ekin, loge);
    }
    else {
      xs = coeff[Z]*ggXsection->GetInelasticElementCrossSection(neutron, ekin,
								Z, aeff[Z]);
    }
  }

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "G4NeutronInelasticXS::ElementCrossSection Z= " << Z
            << " Ekin(MeV)= " << ekin/CLHEP::MeV 
            << ", ElmXSinel(b)= " << xs/CLHEP::barn 
            << G4endl;
  }
#endif
  return xs;
}

G4double
G4NeutronInelasticXS::ComputeIsoCrossSection(G4double ekin, G4double loge,
                                             const G4ParticleDefinition*,
                                             G4int Z, G4int A,
                                             const G4Isotope*, const G4Element*,
                                             const G4Material*)
{
  return IsoCrossSection(ekin, loge, Z, A); 
}

G4double
G4NeutronInelasticXS::GetIsoCrossSection(const G4DynamicParticle* aParticle, 
                                         G4int Z, G4int A,
                                         const G4Isotope*, const G4Element*,
                                         const G4Material*)
{
  return IsoCrossSection(aParticle->GetKineticEnergy(), 
                         aParticle->GetLogKineticEnergy(), Z, A); 
}

G4double
G4NeutronInelasticXS::IsoCrossSection(G4double ekin, G4double loge, 
                                      G4int ZZ, G4int A)
{
  G4double xs{0.0};
  G4int Z = std::min(ZZ, MAXZINEL-1);
  G4bool done{false};

  // data from the resonance region
  if (fRfilesEnabled) {
    auto pv = GetPhysicsVectorR(Z);
    if (nullptr != pv && ekin < inel_max_r_e[Z]) {
      // use isotope x-section if possible
      if (dataR->GetNumberOfComponents(Z) > 0) {
	auto pviso = dataR->GetComponentDataByID(Z, A);
	if (pviso != nullptr) { 
	  xs = pviso->LogVectorValue(ekin, loge);
	  done = true;
	}
      }
      // isotope data are not available or applicable
      if (!done) {
	xs = pv->LogVectorValue(ekin, loge);
	done = true;
      }
    }
  }
  // data above the resonance region
  if (!done) { 
    auto pv = GetPhysicsVector(Z);
    // use isotope x-section if possible
    if (data->GetNumberOfComponents(Z) > 0) {
      auto pviso = data->GetComponentDataByID(Z, A);
      if (pviso != nullptr && ekin <= pviso->GetMaxEnergy()) { 
	xs = pviso->LogVectorValue(ekin, loge);
	done = true;
      }
    }
    // isotope data are not available or applicable
    // use element x-section
    if (!done) {
      if (ekin <= pv->GetMaxEnergy()) {
	xs = pv->LogVectorValue(ekin, loge);
      }
      else {
	xs = coeff[Z]*
	  ggXsection->GetInelasticElementCrossSection(neutron, ekin,
        					      Z, aeff[Z]);
      }
    }
  }

#ifdef G4VERBOSE
  if (verboseLevel > 1) {
    G4cout  << "G4NeutronInelasticXS::IsoXS: Z= " << Z << " A= " << A 
            << " Ekin(MeV)= " << ekin/CLHEP::MeV 
            << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
  return xs;
}

const G4Isotope* G4NeutronInelasticXS::SelectIsotope(
      const G4Element* anElement, G4double kinEnergy, G4double logE)
{
  std::size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);
  if(1 == nIso) { return iso; }

  // more than 1 isotope
  G4int Z = anElement->GetZasInt();
  if (nullptr == data->GetElementData(Z)) { InitialiseOnFly(Z); }

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;
  std::size_t j;

  // isotope wise cross section not available
  if (Z >= MAXZINEL || 0 == data->GetNumberOfComponents(Z)) {
    for (j=0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
        iso = anElement->GetIsotope((G4int)j);
        break;
      }
    }
    return iso;
  }

  // use isotope cross sections
  auto nn = temp.size();
  if(nn < nIso) { temp.resize(nIso, 0.); }

  for (j=0; j<nIso; ++j) {
    sum += abundVector[j]*IsoCrossSection(kinEnergy, logE, Z, 
           anElement->GetIsotope((G4int)j)->GetN());
    temp[j] = sum;
  }
  sum *= q;
  for (j = 0; j<nIso; ++j) {
    if (temp[j] >= sum) {
      iso = anElement->GetIsotope((G4int)j);
      break;
    }
  }
  return iso;
}

void 
G4NeutronInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if (verboseLevel > 0) {
    G4cout << "G4NeutronInelasticXS::BuildPhysicsTable for " 
           << p.GetParticleName() << G4endl;
  }
  if (p.GetParticleName() != "neutron") { 
    G4ExceptionDescription ed;
    ed << p.GetParticleName() << " is a wrong particle type -"
       << " only neutron is allowed";
    G4Exception("G4NeutronInelasticXS::BuildPhysicsTable(..)","had012",
                FatalException, ed, "");
    return; 
  }

  fRfilesEnabled = G4HadronicParameters::Instance()->UseRFilesForXS();

  // it is possible re-initialisation for the new run
  const G4ElementTable* table = G4Element::GetElementTable();

  // initialise static tables only once
  std::call_once(applyOnce, [this]() { isInitializer = true; });

  if (isInitializer) {
    G4AutoLock l(&nInelasticXSMutex);

    // Upload data for elements used in geometry
    for ( auto const & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZINEL-1) );
      if ( nullptr == data->GetElementData(Z) ) { Initialise(Z); }
    }
    l.unlock();
  }
  // prepare isotope selection
  std::size_t nIso = temp.size();
  for ( auto const & elm : *table ) {
    std::size_t n = elm->GetNumberOfIsotopes();
    if (n > nIso) { nIso = n; }
  }
  temp.resize(nIso, 0.0);
}

const G4String& G4NeutronInelasticXS::FindDirectoryPath()
{
  // build the complete string identifying the file with the data set
  if (gDataDirectory.empty()) {
    std::ostringstream ost;
    ost << G4HadronicParameters::Instance()->GetDirPARTICLEXS() << "/neutron/";
    gDataDirectory = ost.str();
  }
  return gDataDirectory;
}

void G4NeutronInelasticXS::InitialiseOnFly(G4int Z)
{
  G4AutoLock l(&nInelasticXSMutex);
  Initialise(Z);
  l.unlock();
}

void G4NeutronInelasticXS::Initialise(G4int Z)
{
  if (nullptr != data->GetElementData(Z)) { return; }

  // upload element data 
  std::ostringstream ost;
  ost << FindDirectoryPath() << "inel" << Z;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);

  G4PhysicsVector* vr = nullptr;
  if (fRfilesEnabled) {
    std::ostringstream ostr;
    ostr << FindDirectoryPath() << "Rinel" << Z;
    vr = RetrieveVector(ostr, false);
    dataR->InitialiseForElement(Z, vr);
  }

  // upload isotope data
  G4bool noComp = true;
  G4bool noCompR = true;
  if (amin[Z] < amax[Z]) {

    for (G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << gDataDirectory << "inel" << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      if (nullptr != v1) {
        if (noComp) {
          G4int nmax = amax[Z] - A + 1;
          data->InitialiseForComponent(Z, nmax);
          noComp = false;
        }
        data->AddComponent(Z, A, v1);
      }
      if (nullptr != vr) {
	std::ostringstream ost2;
	ost2 << gDataDirectory << "Rinel" << Z << "_" << A;
	G4PhysicsVector* v2 = RetrieveVector(ost2, false);
	if (nullptr != v2) {
	  if (noCompR) {
	    G4int nmax = amax[Z] - A + 1;
	    dataR->InitialiseForComponent(Z, nmax);
	    noCompR = false;
	  }
	  dataR->AddComponent(Z, A, v2);
	}
      }
    }
  }
  // no components case
  if (noComp) { data->InitialiseForComponent(Z, 0); }
  if (noCompR) { dataR->InitialiseForComponent(Z, 0); }

  // smooth transition 
  G4double sig1 = (*v)[v->GetVectorLength()-1];
  G4double ehigh= v->GetMaxEnergy();
  G4double sig2 = ggXsection->GetInelasticElementCrossSection(neutron,
                              ehigh, Z, aeff[Z]);
  coeff[Z] = (sig2 > 0.) ? sig1/sig2 : 1.0;
}

G4PhysicsVector* 
G4NeutronInelasticXS::RetrieveVector(std::ostringstream& ost, G4bool warn)
{
  G4PhysicsLogVector* v = nullptr;
  std::ifstream filein(ost.str().c_str());
  if (!filein.is_open()) {
    if (warn) { 
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not opened!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had014",
                  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  } else {
    if (verboseLevel > 1) {
      G4cout << "File " << ost.str() 
             << " is opened by G4NeutronInelasticXS" << G4endl;
    }
    // retrieve data from DB
    v = new G4PhysicsLogVector();
    if (!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
         << "> is not retrieved!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had015",
                  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  return v;
}
