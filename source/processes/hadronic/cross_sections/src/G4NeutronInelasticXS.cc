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
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4IsotopeList.hh"
#include "G4AutoLock.hh"

#include <fstream>
#include <sstream>

G4double G4NeutronInelasticXS::coeff[] = {1.0};
G4ElementData* G4NeutronInelasticXS::data = nullptr;
G4String G4NeutronInelasticXS::gDataDirectory = "";

namespace
{
  G4Mutex nInelasticXSMutex = G4MUTEX_INITIALIZER;
}

G4NeutronInelasticXS::G4NeutronInelasticXS() 
  : G4VCrossSectionDataSet(Default_Name()),
    neutron(G4Neutron::Neutron()),
    elimit(20*CLHEP::MeV)
{
  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout << "G4NeutronInelasticXS::G4NeutronInelasticXS Initialise for Z < " 
	    << MAXZINEL << G4endl;
  }
  ggXsection = G4CrossSectionDataSetRegistry::Instance()->GetComponentCrossSection("Glauber-Gribov");
  if(ggXsection == nullptr) ggXsection = new G4ComponentGGHadronNucleusXsc();
  SetForAllAtomsAndEnergies(true);
}

G4NeutronInelasticXS::~G4NeutronInelasticXS()
{
  if(isMaster) { delete data; data = nullptr; }
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
  return ElementCrossSection(aParticle->GetKineticEnergy(), aParticle->GetLogKineticEnergy(), Z);
}

G4double
G4NeutronInelasticXS::ComputeCrossSectionPerElement(G4double ekin, G4double loge,
						    const G4ParticleDefinition*,
						    const G4Element* elm,
						    const G4Material*)
{
  return ElementCrossSection(ekin, loge, elm->GetZasInt());
}

G4double G4NeutronInelasticXS::ElementCrossSection(G4double ekin, G4double loge, G4int ZZ)
{
  G4int Z = (ZZ >= MAXZINEL) ? MAXZINEL - 1 : ZZ; 
  auto pv = GetPhysicsVector(Z);

  G4double xs = (ekin <= pv->GetMaxEnergy()) ? pv->LogVectorValue(ekin, loge) 
    : coeff[Z]*ggXsection->GetInelasticElementCrossSection(neutron, ekin,
                                                           Z, aeff[Z]);

#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "Z= " << Z << " Ekin(MeV)= " << ekin/CLHEP::MeV 
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
G4NeutronInelasticXS::IsoCrossSection(G4double ekin, G4double logekin, 
                                      G4int ZZ, G4int A)
{
  G4double xs = 0.0;
  G4int Z = (ZZ >= MAXZINEL) ? MAXZINEL - 1 : ZZ;

  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  auto pv = GetPhysicsVector(Z);

  // compute isotope cross section if applicable
  if(ekin <= elimit && amin[Z] < amax[Z] && A >= amin[Z] && A <= amax[Z]) {
    auto pviso = data->GetComponentDataByIndex(Z, A - amin[Z]);
    if(nullptr != pviso) { 
      xs = pviso->LogVectorValue(ekin, logekin); 
#ifdef G4VERBOSE
      if(verboseLevel > 1) {
	G4cout << "G4NeutronInelasticXS::IsoXS: Ekin(MeV)= " 
               << ekin/CLHEP::MeV 
	       << "  xs(b)= " << xs/CLHEP::barn 
	       << "  Z= " << Z << "  A= " << A << G4endl;
      }
#endif
      return xs;
    }
  }
 
  // use element x-section 
  xs = (ekin <= pv->GetMaxEnergy()) ? pv->LogVectorValue(ekin, logekin) 
    : coeff[Z]*ggXsection->GetInelasticElementCrossSection(neutron, ekin,
                                                           Z, aeff[Z]);
  xs *= A/aeff[Z];
#ifdef G4VERBOSE
  if(verboseLevel > 1) {
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
  G4int nIso = (G4int)anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  // more than 1 isotope
  G4int Z = anElement->GetZasInt();
  //G4cout << "SelectIsotope Z= " << Z << G4endl;

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;
  G4int j;

  // isotope wise cross section not available
  if(amax[Z] == amin[Z] || Z >= MAXZINEL) {
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
  G4int nn = (G4int)temp.size();
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
G4NeutronInelasticXS::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if(verboseLevel > 0) {
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

  if(nullptr == data) { 
    G4AutoLock l(&nInelasticXSMutex);
    if(nullptr == data) { 
      isMaster = true;
      data = new G4ElementData();
      data->SetName("NeutronInelastic");
      FindDirectoryPath();
    }
    l.unlock();
  }

  // it is possible re-initialisation for the new run
  const G4ElementTable* table = G4Element::GetElementTable();
  if(isMaster) {
    // Upload data for elements used in geometry
    for ( auto & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZINEL-1) );
      if ( nullptr == data->GetElementData(Z) ) { Initialise(Z); }
    }
  }
  // prepare isotope selection
  std::size_t nIso = temp.size();
  for ( auto & elm : *table ) {
    std::size_t n = elm->GetNumberOfIsotopes();
    if(n > nIso) { nIso = n; }
  }
  temp.resize(nIso, 0.0);
}

const G4String& G4NeutronInelasticXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4PARTICLEXSDATA");
    if (nullptr != path) {
      std::ostringstream ost;
      ost << path << "/neutron/inel";
      gDataDirectory = ost.str();
    } else {
      G4Exception("G4NeutronInelasticXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4PARTICLEXSDATA is not defined");
    }
  }
  return gDataDirectory;
}

void G4NeutronInelasticXS::InitialiseOnFly(G4int Z)
{
  if(nullptr == data->GetElementData(Z)) {
    G4AutoLock l(&nInelasticXSMutex);
    if(nullptr == data->GetElementData(Z)) { Initialise(Z); }
    l.unlock();
  }
}

void G4NeutronInelasticXS::Initialise(G4int Z)
{
  if(nullptr != data->GetElementData(Z)) { return; }

  // upload element data 
  std::ostringstream ost;
  ost <<  FindDirectoryPath() << Z;
  G4PhysicsVector* v = RetrieveVector(ost, true);
  data->InitialiseForElement(Z, v);
  if(verboseLevel > 1) {
    G4cout  << "G4NeutronInelasticXS::Initialise for Z= " << Z 
	    << " A= " << aeff[Z] << "  Amin= " << amin[Z] 
	    << "  Amax= " << amax[Z] << G4endl;
  }
  // upload isotope data
  if(amin[Z] < amax[Z]) {
    G4int nmax = amax[Z] - amin[Z] + 1;
    data->InitialiseForComponent(Z, nmax);

    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << gDataDirectory << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false);
      data->AddComponent(Z, A, v1); 
    }
  }

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
    if(warn) { 
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4NeutronInelasticXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
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
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
  return v;
}
