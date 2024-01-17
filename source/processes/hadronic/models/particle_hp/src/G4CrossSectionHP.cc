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
// V. Ivanchenko September 2023 
//               
// G4CrossSectionHP is a generic class implementing 
// cross sections for neutrons, protons and light ions
// It is an alternative to code developed by J.P. Wellisch & T.Koi
//

#include <fstream>
#include <sstream>
#include <thread>

#include "G4CrossSectionHP.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementDataRegistry.hh"
#include "G4ElementTable.hh"
#include "G4IsotopeList.hh"
#include "G4HadronicParameters.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex theHPXS = G4MUTEX_INITIALIZER;
}

G4CrossSectionHP::G4CrossSectionHP(const G4ParticleDefinition* p,
                                   const G4String& nameData,
                                   const G4String& nameDir, G4double emaxHP,
                                   G4int zmin, G4int zmax)
  : G4VCrossSectionDataSet(nameData),
    fParticle(p),
    fManagerHP(G4ParticleHPManager::GetInstance()),
    emax(emaxHP),
    emaxT(fManagerHP->GetMaxEnergyDoppler()),
    elimit(1.0e-11*CLHEP::eV),
    logElimit(G4Log(elimit)),
    minZ(zmin),
    maxZ(zmax),
    fDataName(nameData),
    fDataDirectory(nameDir)
{
  //  verboseLevel = 1;
  if (verboseLevel > 1) {
    G4cout << "G4CrossSectionHP::G4CrossSectionHP: Initialise for "
	   << fDataName << "  " << minZ << " < Z < " << maxZ
	   << "  EmaxT(MeV)=" << emaxT << G4endl;
    G4cout << "Data directory: " << fDataDirectory << G4endl; 
  }
  auto ptr = G4ElementDataRegistry::Instance();
  auto data = ptr->GetElementDataByName(fDataName);
  if (nullptr == data) { 
    data = new G4ElementData(maxZ - minZ + 1);
    data->SetName(fDataName);
  }
  fData = data;
}

G4bool G4CrossSectionHP::IsIsoApplicable(const G4DynamicParticle* dp,
					 G4int, G4int,
                                         const G4Element*,
					 const G4Material*)
{
  return (dp->GetKineticEnergy() <= emax);
}

G4double G4CrossSectionHP::GetIsoCrossSection(const G4DynamicParticle* dp,
                                              G4int Z, G4int A,
                                              const G4Isotope*,
					      const G4Element*,
					      const G4Material* mat)
{
  G4double ekin = dp->GetKineticEnergy();
  G4double loge = dp->GetLogKineticEnergy();
  if (ekin < elimit) {
    ekin = elimit;
    loge = logElimit;
  }
  if (mat != fCurrentMat) { PrepareCache(mat); }

  return IsoCrossSection(ekin, loge, Z, A, mat->GetTemperature());
}

G4double
G4CrossSectionHP::ComputeIsoCrossSection(G4double e, G4double le,
                                         const G4ParticleDefinition*,
                                         G4int Z, G4int A,
					 const G4Isotope*,
					 const G4Element*,
					 const G4Material* mat)
{
  G4double ekin = e;
  G4double loge = le;
  if (ekin < elimit) {
    ekin = elimit;
    loge = logElimit;
  }
  if (mat != fCurrentMat) { PrepareCache(mat); }

  return IsoCrossSection(ekin, loge, Z, A, mat->GetTemperature()); 
}

G4double G4CrossSectionHP::IsoCrossSection(const G4double ekin,
                                           const G4double logek,
					   const G4int Z, const G4int A,
                                           const G4double T)
{
  //G4cout << "G4CrossSectionHP::IsoCrossSection Z=" << Z << " A=" << A
  // << " E(MeV)=" << ekin/MeV << " T=" << T << "  " << GetName() << G4endl;
  G4double xs = 0.0;
  if (ekin > emax || Z > maxZ || Z < minZ || ekin < elimit) { return xs; }

  const G4PhysicsVector* pv0 = fData->GetElementData(Z - minZ);
  if (nullptr == pv0) {
    InitialiseOnFly(Z);
    pv0 = fData->GetElementData(Z - minZ);
  }
  if (nullptr == pv0) { return xs; }
  const G4PhysicsVector* pv = fData->GetComponentDataByID(Z - minZ, A);
  if (nullptr == pv) { return xs; }

  // no Doppler broading
  G4double factT = T/CLHEP::STP_Temperature;
  if (ekin >= emaxT*factT || fManagerHP->GetNeglectDoppler()) {
    xs = pv->LogFreeVectorValue(ekin, logek);

  } else {

    // Doppler broading
    G4double lambda = 1.0/(CLHEP::k_Boltzmann*T);
    G4double mass = fParticle->GetPDGMass();
    G4double massTarget = G4NucleiProperties::GetNuclearMass(A, Z);
    G4LorentzVector lv(0., 0., 0., mass + ekin);

    // limits of integration
    const G4double lim = 1.01;
    const G4int nmin = 3;
    G4int i;
    const G4int nn = 20;
    G4double xs2 = 0.0;

    for (i=1; i<nn; ++i) {
      G4double erand = G4RandGamma::shoot(2.0, lambda);
      auto mom = G4RandomDirection()*std::sqrt(2*massTarget*erand);
      fLV.set(mom.x(), mom.y(), mom.z(), mass + erand);
      fBoost = fLV.boostVector();
      G4double e = lv.boost(fBoost).e() - mass;
      G4double y = pv->Value(e, index);
      xs += y;
      xs2 += y*y;
      if (i >= nmin && i*xs2 <= lim*xs*xs) { break; } 
    }
    xs /= (G4double)std::min(i, nn-1);
  }
#ifdef G4VERBOSE
  if (verboseLevel > 1) {
    G4cout << "G4CrossSectionHP::IsoXS " << fDataName 
	   << "  Z= " << Z << "  A= " << A 
	   << " Ekin(MeV)= " << ekin/MeV << "  xs(b)= " << xs/barn 
           << "  size=" << fZA.size() << G4endl;
  }
#endif

  // save cross section into struct
  for (std::size_t i=0; i<fZA.size(); ++i) {
    if (Z == fZA[i].first && A == fZA[i].second) {
      fIsoXS[i] = xs;
      break;
    }
  }
  return xs;
}

const G4Isotope* G4CrossSectionHP::SelectIsotope(const G4Element* elm,
                                                 G4double, G4double)
{
  std::size_t nIso = elm->GetNumberOfIsotopes();
  const G4Isotope* iso = elm->GetIsotope(0);

  //G4cout << "SelectIsotope NIso= " << nIso << G4endl;
  if(1 == nIso) { return iso; }

  // more than 1 isotope
  G4int Z = elm->GetZasInt();
  if (Z >= minZ && Z <= maxZ && nullptr == fData->GetElementData(Z - minZ)) {
    InitialiseOnFly(Z);
  }
  
  const G4double* abundVector = elm->GetRelativeAbundanceVector();
  G4double q = G4UniformRand();
  G4double sum = 0.0;

  // is there cached isotope wise cross section?
  std::size_t j;
  if (Z < minZ || Z > maxZ || !CheckCache(Z) ||
      0 == fData->GetNumberOfComponents(Z - minZ)) {
    for (j = 0; j<nIso; ++j) {
      sum += abundVector[j];
      if(q <= sum) {
	iso = elm->GetIsotope((G4int)j);
	break;
      }
    }
    return iso;
  } 
  std::size_t nn = fTemp.size();
  if (nn < nIso) { fTemp.resize(nIso, 0.); }

  // reuse cache
  for (j=0; j<nIso; ++j) {
    sum += abundVector[j]*
      GetCrossSection(Z - minZ, elm->GetIsotope((G4int)j)->GetN());
    fTemp[j] = sum;
  }
  sum *= q;
  for (j = 0; j<nIso; ++j) {
    if (fTemp[j] >= sum) {
      iso = elm->GetIsotope((G4int)j);
      break;
    }
  }
  return iso;
}

void G4CrossSectionHP::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  if (verboseLevel > 1){
    G4cout << "G4CrossSectionHP::BuildPhysicsTable for "
	   << p.GetParticleName() << " and " << fDataName << G4endl;
  }

  // it is possible re-initialisation for the second run
  const G4ElementTable* table = G4Element::GetElementTable();

  // Access to elements
  for ( auto const & elm : *table ) {
    G4int Z = elm->GetZasInt();
    if (Z >= minZ && Z <= maxZ &&
        nullptr == fData->GetElementData(Z - minZ)) { 
      InitialiseOnFly(Z);
    }
  }

  // prepare isotope selection
  std::size_t nmax = 0;
  std::size_t imax = 0;
  for ( auto const & mat : *G4Material::GetMaterialTable() ) {
    std::size_t n = 0; 
    for ( auto const & elm : *mat->GetElementVector() ) {
      std::size_t niso = elm->GetNumberOfIsotopes();
      n += niso;
      if(niso > imax) { imax = niso; }
    }
    if (n > nmax) { nmax = n; }
  }
  fTemp.resize(imax, 0.0);
  fZA.clear();
  fZA.reserve(nmax);
  fIsoXS.resize(nmax, 0.0);
}

void G4CrossSectionHP::DumpPhysicsTable(const G4ParticleDefinition&)
{
  if (fManagerHP->GetVerboseLevel() == 0 || fPrinted)
    return;

  //
  // Dump element based cross section
  // range 10e-5 eV to 20 MeV
  // 10 point per decade
  // in barn
  //
  fPrinted = true;
  G4cout << G4endl;
  G4cout << "HP Cross Section " << fDataName << " for "
	 << fParticle->GetParticleName() << G4endl;
  G4cout << "(Pointwise cross-section at 0 Kelvin.)" << G4endl;
  G4cout << G4endl;
  G4cout << "Name of Element" << G4endl;
  G4cout << "Energy[eV]  XS[barn]" << G4endl;
  G4cout << G4endl;

  const G4ElementTable* table = G4Element::GetElementTable();
  for ( auto const & elm : *table ) {
    G4int Z = elm->GetZasInt();
    if (Z >= minZ && Z <= maxZ &&
        nullptr != fData->GetElementData(Z - minZ)) {
      G4cout << "---------------------------------------------------" << G4endl;
      G4cout << elm->GetName() << G4endl;
      std::size_t n = fData->GetNumberOfComponents(Z);
      for (size_t i=0; i<n; ++i) {
	G4cout << "##  Z=" << Z << "  A=" << fData->GetComponentID(Z - minZ, i);
	G4cout << *(fData->GetComponentDataByIndex(Z - minZ, i)) << G4endl;
      }
    }
  }
}

void G4CrossSectionHP::PrepareCache(const G4Material* mat)
{
  fCurrentMat = mat;
  fZA.clear();
  for ( auto const & elm : *(mat->GetElementVector()) ) {
    G4int Z = elm->GetZasInt();
    for ( auto const & iso : *(elm->GetIsotopeVector()) ) {
      G4int A = iso->GetN();
      fZA.emplace_back(Z, A);
    }
  }
  fIsoXS.resize(fZA.size(), 0.0);
}

void G4CrossSectionHP::InitialiseOnFly(const G4int Z)
{
  G4AutoLock l(&theHPXS);
  Initialise(Z);
  l.unlock();
}

void G4CrossSectionHP::Initialise(const G4int Z)
{
  if (fManagerHP->GetVerboseLevel() > 1) {
    G4cout << " G4CrossSectionHP::Initialise: Z=" << Z 
           << " for " << fDataName
	   << " minZ=" << minZ << " maxZ=" << maxZ << G4endl;
  }
  if (Z < minZ || Z > maxZ || nullptr != fData->GetElementData(Z - minZ)) { 
    return;
  }

  // add empty vector to avoid double initialisation
  fData->InitialiseForElement(Z - minZ, new G4PhysicsVector());

  G4String tnam = "temp";
  G4bool noComp = true;
  for (G4int A=amin[Z]; A<=amax[Z]; ++A) {
    std::ostringstream ost;
    ost << fDataDirectory << Z << "_";
    if (6 == Z && 12 == A) {
      ost << "nat_";
    } else if (27 == Z && 62 == A) {
      ost << "62m1_";
    } else if (47 == Z && 106 == A) {
      ost << "106m1_";
    } else if (48 == Z && 115 == A) {
      ost << "115m1_";
    } else if (52 == Z && 127 == A) {
      ost << "127m1_";
    } else if (52 == Z && 129 == A) {
      ost << "129m1_";
    } else if (52 == Z && 131 == A) {
      ost << "131m1_";
    } else if (67 == Z && 166 == A) {
      ost << "166m1_";
    } else if (73 == Z && 180 == A) {
      ost << "180m1_";
    } else {
      ost << A << "_";
    }
    ost << elementName[Z];
    std::ifstream filein(ost.str().c_str());
    //G4cout << "File: " << ost.str() << "  " << G4endl;
    std::istringstream theXSData(tnam, std::ios::in);
    fManagerHP->GetDataStream(ost.str().c_str(), theXSData);
    if (theXSData) {
      G4int i1, i2, n;
      theXSData >> i1 >> i2 >> n;
      if (fManagerHP->GetVerboseLevel() > 1) {
	G4cout << "## G4CrossSectionHP::Initialise for Z=" << Z
	       << " A=" << A << " Npoints=" << n << G4endl;
      }
      G4double x, y;
      G4PhysicsFreeVector* v = new G4PhysicsFreeVector(n);
      for (G4int i=0; i<n; ++i) {
	theXSData >> x >> y;
	x *= CLHEP::eV;
	y *= CLHEP::barn;         
	//G4cout << "  e=" << x << "  xs=" << y << G4endl;
	v->PutValues((std::size_t)i, x, y);
      }
      v->EnableLogBinSearch(2);
      if (noComp) {
	G4int nmax = amax[Z] - A + 1;
	fData->InitialiseForComponent(Z - minZ, nmax);
	noComp = false;
      }
      fData->AddComponent(Z - minZ, A, v);
    }
  }   
  if (noComp) { fData->InitialiseForComponent(Z - minZ, 0); }
}

