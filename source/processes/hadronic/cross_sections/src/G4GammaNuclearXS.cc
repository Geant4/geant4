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
// Authors  V.Ivantchenko, Geant4, 20 October 2020
//          B.Kutsenko, BINP/NSU, 10 August 2021
//
// Modifications:
//

#include "G4GammaNuclearXS.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4Gamma.hh"
#include "G4IsotopeList.hh"

#include <fstream>
#include <sstream>
#include <vector>

G4ElementData* G4GammaNuclearXS::data = nullptr;

G4double G4GammaNuclearXS::coeff[3][3];
G4double G4GammaNuclearXS::xs150[MAXZGAMMAXS] = {0.0};
G4String G4GammaNuclearXS::gDataDirectory = "";

#ifdef G4MULTITHREADED
  G4Mutex G4GammaNuclearXS::gNuclearXSMutex = G4MUTEX_INITIALIZER;
#endif

G4GammaNuclearXS::G4GammaNuclearXS() 
  : G4VCrossSectionDataSet(Default_Name()),
   gamma(G4Gamma::Gamma())
{
  //  verboseLevel = 0;
  if(verboseLevel > 0){
    G4cout  << "G4GammaNuclearXS::G4GammaNuclearXS Initialise for Z < " 
	    << MAXZGAMMAXS << G4endl;
  }
  ggXsection = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet("PhotoNuclearXS");
  if(ggXsection == nullptr) ggXsection = new G4PhotoNuclearCrossSection();
  SetForAllAtomsAndEnergies(true);
}

G4GammaNuclearXS::~G4GammaNuclearXS()
{
  if(isMaster) { 
    delete data; 
    data = nullptr; 
  }
}

void G4GammaNuclearXS::CrossSectionDescription(std::ostream& outFile) const
{
  outFile << "G4GammaNuclearXS calculates the gamma nuclear\n"
          << "cross-section for GDR energy region on nuclei using data from the high precision\n"
          << "IAEA photonuclear database (2019). Then liniear connection\n"
	  <<"implemented with previous CHIPS photonuclear model\n";
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
  const G4int Z = (ZZ >= MAXZGAMMAXS) ? MAXZGAMMAXS - 1 : ZZ;
  auto pv = GetPhysicsVector(Z);

  if(pv == nullptr) {
    return ggXsection->GetElementCrossSection(aParticle, Z, mat);
  }
  const G4double emax = pv->GetMaxEnergy();
  const G4double ekin = aParticle->GetKineticEnergy();
  G4double xs = 0.0;
  if(ekin <= emax) {
    xs = pv->Value(ekin);
  } else if(ekin >= rTransitionBound){
    xs = ggXsection->GetElementCrossSection(aParticle, Z, mat);
  } else {
    const G4double rxs = xs150[Z];
    const G4double lxs = pv->Value(emax);
    xs = lxs + (ekin - emax)*(rxs - lxs)/(rTransitionBound-emax);
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

G4double 
G4GammaNuclearXS::ElementCrossSection(G4double ekin, G4int ZZ)
{    
  G4DynamicParticle theGamma(gamma, G4ThreeVector(1,0,0), ekin);
  return GetElementCrossSection(&theGamma, ZZ);
}

G4double 
G4GammaNuclearXS::IsoCrossSection(G4double ekin, G4int ZZ, G4int A)
{
  G4DynamicParticle theGamma(gamma, G4ThreeVector(1,0,0), ekin);
  return GetIsoCrossSection(&theGamma, ZZ, A);
}

G4double G4GammaNuclearXS::GetIsoCrossSection(
         const G4DynamicParticle* aParticle,
	 G4int ZZ, G4int A,
	 const G4Isotope*, const G4Element*,
         const G4Material*)
{
  const G4int Z = (ZZ >= MAXZGAMMAXS) ? MAXZGAMMAXS - 1 : ZZ; 
  /*
  G4cout << "IsoCrossSection  Z= " << Z << "  A= " << A 
         << "  Amin= " << amin[Z] << " Amax= " << amax[Z]
         << " E(MeV)= " << ekin << G4endl;
  */
  auto pv = GetPhysicsVector(Z);
  if(pv == nullptr) {
    return ggXsection->GetIsoCrossSection(aParticle, Z, A);
  }
  const G4double ekin = aParticle->GetKineticEnergy();
  const G4double emax = pv->GetMaxEnergy();
  G4double xs = 0.0;

  // compute isotope cross section if applicable
  if(amin[Z] < amax[Z] && A >= amin[Z] && A <= amax[Z] && 
     ekin < rTransitionBound) {
    auto pviso = data->GetComponentDataByIndex(Z, A - amin[Z]);
    // isotope file exists
    if(nullptr != pviso) {
      const G4double emaxiso = pviso->GetMaxEnergy();
      if(ekin <= emaxiso) {
	xs = pviso->Value(ekin);
      } else {
	G4DynamicParticle 
	  theGamma(gamma, G4ThreeVector(0,0,1.), rTransitionBound);
	const G4double rxs = ggXsection->GetIsoCrossSection(&theGamma, Z, A);
	const G4double lxs = pviso->Value(emaxiso);
	xs = lxs + (ekin - emaxiso)*(rxs - lxs)/(rTransitionBound-emaxiso);
      }   
#ifdef G4VERBOSE
      if(verboseLevel > 1) {
	G4cout  << "G4GammaNuclearXS::IsoXS: Z= " << Z << " A= " << A 
		<< " Ekin(MeV)= " << ekin/CLHEP::MeV 
		<< ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
      }
#endif
      return xs;
    }
  }

  // use element x-section
  // for the hydrogen target there is no element data
  if(ekin <= emax && Z != 1) { 
    xs = pv->Value(ekin)*A/aeff[Z];

    // CHIPS for high energy and for the hydrogen target
  } else if(ekin >= rTransitionBound || Z == 1) {
    if(Z <= 2 && ekin > 10.*GeV) { 
      xs = coeff[Z][A - amin[Z]]*
	ggXsection->GetElementCrossSection(aParticle, Z, 0);
    } else {
      xs = ggXsection->GetIsoCrossSection(aParticle, Z, A);
    }

    // transition GDR to CHIPS
  } else {
    const G4double rxs = xs150[Z];
    const G4double lxs = pv->Value(emax)*A/aeff[Z];
    xs = lxs + (ekin - emax)*(rxs - lxs)/(rTransitionBound-emax);
  }
#ifdef G4VERBOSE
  if(verboseLevel > 1) {
    G4cout  << "G4GammaNuclearXS::IsoXS: Z= " << Z << " A= " << A 
	    << " Ekin(MeV)= " << ekin/CLHEP::MeV 
	    << ", ElmXS(b)= " << xs/CLHEP::barn << G4endl;
  }
#endif
  return xs;
}

const G4Isotope* G4GammaNuclearXS::SelectIsotope(
       const G4Element* anElement, G4double kinEnergy, G4double)
{
  std::size_t nIso = anElement->GetNumberOfIsotopes();
  const G4Isotope* iso = anElement->GetIsotope(0);

  if(1 == nIso) { return iso; }

  const G4double* abundVector = anElement->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;
  G4int j;
  G4int Z = anElement->GetZasInt();

  // condition to use only isotope abundance
  if(amax[Z] == amin[Z] || kinEnergy > rTransitionBound || Z >= MAXZGAMMAXS ) {
    for (j=0; j<(G4int)nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = anElement->GetIsotope(j);
	break;
      }
    }
    return iso;
  }
  // use isotope cross sections
  std::size_t nn = temp.size();
  if(nn < nIso) { temp.resize(nIso, 0.); }
  
  for (j=0; j<(G4int)nIso; ++j) {
    //G4cout << j << "-th isotope " << (*isoVector)[j]->GetN() 
    //       <<  " abund= " << abundVector[j] << G4endl;
    sum += abundVector[j]*
      IsoCrossSection(kinEnergy, Z, anElement->GetIsotope(j)->GetN());
    temp[j] = sum;
  }
  sum *= q;
  for (j = 0; j<(G4int)nIso; ++j) {
    if(temp[j] >= sum) {
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

    if(nullptr == data) { 
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&gNuclearXSMutex);
    if(nullptr == data) { 
#endif
      isMaster = true;
      data = new G4ElementData(); 
      data->SetName("PhotoNuclear");
      FindDirectoryPath();
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&gNuclearXSMutex);
#endif
  }
    

  // it is possible re-initialisation for the second run
  // Upload data for elements used in geometry
  const G4ElementTable* table = G4Element::GetElementTable();
  if(isMaster) {
    for ( auto & elm : *table ) {
      G4int Z = std::max( 1, std::min( elm->GetZasInt(), MAXZGAMMAXS-1) );
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

const G4String& G4GammaNuclearXS::FindDirectoryPath()
{
  // check environment variable
  // build the complete string identifying the file with the data set
  if(gDataDirectory.empty()) {
    const char* path = G4FindDataDir("G4PARTICLEXSDATA");
    if (nullptr != path) {
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
   if(nullptr == data->GetElementData(Z)) { 
#endif
     Initialise(Z);
#ifdef G4MULTITHREADED
   }
   G4MUTEXUNLOCK(&gNuclearXSMutex);
#endif
}

void G4GammaNuclearXS::Initialise(G4int Z)
{
  if(nullptr != data->GetElementData(Z)) { return; }

  // upload data from file
  std::ostringstream ost;
  ost << FindDirectoryPath() << Z ;
  G4PhysicsVector* v = RetrieveVector(ost, true, Z);
  
  data->InitialiseForElement(Z, v);
  /*
  G4cout << "G4GammaNuclearXS::Initialise for Z= " << Z 
	 << " A= " << Amean << "  Amin= " << amin[Z] 
	 << "  Amax= " << amax[Z] << G4endl;
  */
  // upload isotope data
  G4DynamicParticle theGamma(gamma, G4ThreeVector(1,0,0), rTransitionBound);
  xs150[Z] = ggXsection->GetElementCrossSection(&theGamma, Z, 0);
  if(amax[Z] > amin[Z]) {
    G4int nmax = amax[Z]-amin[Z]+1;
    data->InitialiseForComponent(Z, nmax);
    for(G4int A=amin[Z]; A<=amax[Z]; ++A) {
      std::ostringstream ost1;
      ost1 << gDataDirectory << Z << "_" << A;
      G4PhysicsVector* v1 = RetrieveVector(ost1, false, Z);
      data->AddComponent(Z, A, v1);
      if(Z<=2){
	theGamma.SetKineticEnergy(10.*GeV);
	G4double sig1 = ggXsection->GetIsoCrossSection(&theGamma, Z, A);
	G4double sig2 = ggXsection->GetElementCrossSection(&theGamma, Z, 0);
	if(sig2 > 0.) coeff[Z][A-amin[Z]]=(sig1/sig2);
	else coeff[Z][A-amin[Z]]=1.;
      }
    }
  }
}

G4PhysicsVector* 
G4GammaNuclearXS::RetrieveVector(std::ostringstream& ost, G4bool isElement, G4int Z)
{
  G4PhysicsVector* v = nullptr;

  std::ifstream filein(ost.str().c_str());
  if (!filein.is_open()) {
    if(isElement) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not opened!";
      G4Exception("G4GammaNuclearXS::RetrieveVector(..)","had014",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  } else {
    if(verboseLevel > 1) {
      G4cout << "File " << ost.str() 
	     << " is opened by G4GammaNuclearXS" << G4endl;
    }
    // retrieve data from DB
    if(std::find(std::begin(freeVectorException), std::end(freeVectorException), Z ) == std::end(freeVectorException) && isElement) {
      v = new G4PhysicsLinearVector();
    } else {
      v = new G4PhysicsVector();
    }
    if(!v->Retrieve(filein, true)) {
      G4ExceptionDescription ed;
      ed << "Data file <" << ost.str().c_str()
	 << "> is not retrieved!";
      G4Exception("G4GammaNuclearXS::RetrieveVector(..)","had015",
		  FatalException, ed, "Check G4PARTICLEXSDATA");
    }
  }
 return v;
}

