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
// this code implementation is the intellectual property of
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard (and T. Koi)
// 081203 limit maximum trial for creating final states add protection for 1H isotope case by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//
#include "G4ParticleHPInelastic.hh"

#include "G4HadronicParameters.hh"
#include "G4ParticleHP2AInelasticFS.hh"
#include "G4ParticleHP2N2AInelasticFS.hh"
#include "G4ParticleHP2NAInelasticFS.hh"
#include "G4ParticleHP2NDInelasticFS.hh"
#include "G4ParticleHP2NInelasticFS.hh"
#include "G4ParticleHP2NPInelasticFS.hh"
#include "G4ParticleHP2PInelasticFS.hh"
#include "G4ParticleHP3AInelasticFS.hh"
#include "G4ParticleHP3NAInelasticFS.hh"
#include "G4ParticleHP3NInelasticFS.hh"
#include "G4ParticleHP3NPInelasticFS.hh"
#include "G4ParticleHP4NInelasticFS.hh"
#include "G4ParticleHPAInelasticFS.hh"
#include "G4ParticleHPD2AInelasticFS.hh"
#include "G4ParticleHPDAInelasticFS.hh"
#include "G4ParticleHPDInelasticFS.hh"
#include "G4ParticleHPHe3InelasticFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPN2AInelasticFS.hh"
#include "G4ParticleHPN2PInelasticFS.hh"
#include "G4ParticleHPN3AInelasticFS.hh"
#include "G4ParticleHPNAInelasticFS.hh"
#include "G4ParticleHPND2AInelasticFS.hh"
#include "G4ParticleHPNDInelasticFS.hh"
#include "G4ParticleHPNHe3InelasticFS.hh"
#include "G4ParticleHPNInelasticFS.hh"
#include "G4ParticleHPNPAInelasticFS.hh"
#include "G4ParticleHPNPInelasticFS.hh"
#include "G4ParticleHPNT2AInelasticFS.hh"
#include "G4ParticleHPNTInelasticFS.hh"
#include "G4ParticleHPNXInelasticFS.hh"
#include "G4ParticleHPPAInelasticFS.hh"
#include "G4ParticleHPPDInelasticFS.hh"
#include "G4ParticleHPPInelasticFS.hh"
#include "G4ParticleHPPTInelasticFS.hh"
#include "G4ParticleHPT2AInelasticFS.hh"
#include "G4ParticleHPTInelasticFS.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

G4bool G4ParticleHPInelastic::fLock[] = {true, true, true, true, true, true};
std::vector<G4ParticleHPChannelList*>*
G4ParticleHPInelastic::theInelastic[] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

namespace
{
  G4Mutex theHPInelastic = G4MUTEX_INITIALIZER;
}

G4ParticleHPInelastic::G4ParticleHPInelastic(G4ParticleDefinition* p, const char* name)
  : G4HadronicInteraction(name), theProjectile(p)
{
  fManager = G4ParticleHPManager::GetInstance();
  dirName = fManager->GetParticleHPPath(theProjectile) + "/Inelastic";
  indexP = fManager->GetPHPIndex(theProjectile);

#ifdef G4VERBOSE
  if (fManager->GetVerboseLevel() > 1)
    G4cout << "@@@ G4ParticleHPInelastic instantiated for "
           << p->GetParticleName() << " indexP=" << indexP
	   << "/n    data directory " << dirName << G4endl;
#endif
}

G4ParticleHPInelastic::~G4ParticleHPInelastic()
{
  // Vector is shared, only one delete
  if (isFirst) {
    ClearData();
  }
}

void G4ParticleHPInelastic::ClearData()
{
  if (theInelastic[indexP] != nullptr) {
    for (auto const& p : *(theInelastic[indexP])) {
      delete p;
    }
    delete theInelastic[indexP];
    theInelastic[indexP] = nullptr;
  }
}

G4HadFinalState* G4ParticleHPInelastic::ApplyYourself(const G4HadProjectile& aTrack,
                                                      G4Nucleus& aNucleus)
{
  G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
  const G4Material* theMaterial = aTrack.GetMaterial();
  auto n = (G4int)theMaterial->GetNumberOfElements();
  auto elm = theMaterial->GetElement(0);
  std::size_t index = elm->GetIndex();
  G4int it = 0;
  /*
  G4cout << "G4ParticleHPInelastic::ApplyYourself n=" << n << " index=" << index
	 << " indexP=" << indexP << " "
         << aTrack.GetDefinition()->GetParticleName() << G4endl;
  */
  if (n != 1) {
    auto xSec = new G4double[n];
    G4double sum = 0;
    G4int i;
    const G4double* NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;
    G4double xs; 
    G4ParticleHPThermalBoost aThermalE;
    for (i = 0; i < n; ++i) {
      elm = theMaterial->GetElement(i);
      index = elm->GetIndex();
      /*
      G4cout << "i=" << i << "  index=" << index << "  " << elm->GetName() 
	     << "  " << (*(theInelastic[indexP]))[index] << G4endl;
      */
      rWeight = NumAtomsPerVolume[i];
      if (aTrack.GetDefinition() == G4Neutron::Neutron()) {
        xs = (*(theInelastic[indexP]))[index]->GetXsec(aThermalE.GetThermalEnergy(aTrack, elm,
						       theMaterial->GetTemperature()));
      }
      else {
        xs = (*(theInelastic[indexP]))[index]->GetXsec(aTrack.GetKineticEnergy());
      }
      xs *= rWeight;
      sum += xs;
      xSec[i] = sum;
#ifdef G4VERBOSE
      if (fManager->GetDEBUG())
        G4cout << " G4ParticleHPInelastic XSEC ELEM " << i << " = " << xSec[i] << G4endl;
#endif
    }
    sum *= G4UniformRand();
    for (it = 0; it < n; ++it) {
      elm = theMaterial->GetElement(it);
      index = elm->GetIndex();
      if (sum <= xSec[it]) break;
    }
    delete[] xSec;
  }

#ifdef G4VERBOSE
  if (fManager->GetDEBUG())
    G4cout << " G4ParticleHPInelastic: Elem it=" << it << "  "
           << elm->GetName() << " index=" << index
	   << " from material " << theMaterial->GetName()
           << G4endl;
#endif

  G4HadFinalState* result = 
    (*(theInelastic[indexP]))[index]->ApplyYourself(elm, aTrack);

  aNucleus.SetParameters(fManager->GetReactionWhiteBoard()->GetTargA(),
                         fManager->GetReactionWhiteBoard()->GetTargZ());
  
  const G4Element* target_element = (*G4Element::GetElementTable())[index];
  const G4Isotope* target_isotope = nullptr;
  auto iele = (G4int)target_element->GetNumberOfIsotopes();
  for (G4int j = 0; j != iele; ++j) {
    target_isotope = target_element->GetIsotope(j);
    if (target_isotope->GetN() == fManager->GetReactionWhiteBoard()->GetTargA())
      break;
  }
  aNucleus.SetIsotope(target_isotope);

  G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();

  return result;
}

const std::pair<G4double, G4double> G4ParticleHPInelastic::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus
  return std::pair<G4double, G4double>(10.0 * perCent, 350.0 * CLHEP::GeV);
}

void G4ParticleHPInelastic::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (fLock[indexP]) {
    G4AutoLock l(&theHPInelastic);
    if (fLock[indexP]) {
      isFirst = true;
      fLock[indexP] = false;
    }
    l.unlock();
  }

  G4int nelm = (G4int)G4Element::GetNumberOfElements();
  G4int n0 = numEle;
  numEle = nelm;
  if (!isFirst || nelm == n0) { return; }

  // extra elements should be initialized
  G4AutoLock l(&theHPInelastic);

  if (nullptr == theInelastic[indexP]) {
    theInelastic[indexP] = new std::vector<G4ParticleHPChannelList*>;
  }

  if (fManager->GetVerboseLevel() > 0 && isFirst) {
    fManager->DumpSetting();
    G4cout << "@@@ G4ParticleHPInelastic instantiated for particle "
	   << theProjectile->GetParticleName() << "/n    data directory is "
	   << dirName << G4endl;
  }

  auto table = G4Element::GetElementTable();
  for (G4int i = n0; i < nelm; ++i) {
    auto clist = new G4ParticleHPChannelList(36, theProjectile);
    theInelastic[indexP]->push_back(clist);
    clist->Init((*table)[i], dirName, theProjectile);
    clist->Register(new G4ParticleHPNInelasticFS, "F01/");  // has
    clist->Register(new G4ParticleHPNXInelasticFS, "F02/");
    clist->Register(new G4ParticleHP2NDInelasticFS, "F03/");
    clist->Register(new G4ParticleHP2NInelasticFS, "F04/");  // has, E Done
    clist->Register(new G4ParticleHP3NInelasticFS, "F05/");  // has, E Done
    clist->Register(new G4ParticleHPNAInelasticFS, "F06/");
    clist->Register(new G4ParticleHPN3AInelasticFS, "F07/");
    clist->Register(new G4ParticleHP2NAInelasticFS, "F08/");
    clist->Register(new G4ParticleHP3NAInelasticFS, "F09/");
    clist->Register(new G4ParticleHPNPInelasticFS, "F10/");
    clist->Register(new G4ParticleHPN2AInelasticFS, "F11/");
    clist->Register(new G4ParticleHP2N2AInelasticFS, "F12/");
    clist->Register(new G4ParticleHPNDInelasticFS, "F13/");
    clist->Register(new G4ParticleHPNTInelasticFS, "F14/");
    clist->Register(new G4ParticleHPNHe3InelasticFS, "F15/");
    clist->Register(new G4ParticleHPND2AInelasticFS, "F16/");
    clist->Register(new G4ParticleHPNT2AInelasticFS, "F17/");
    clist->Register(new G4ParticleHP4NInelasticFS, "F18/");  // has, E Done
    clist->Register(new G4ParticleHP2NPInelasticFS, "F19/");
    clist->Register(new G4ParticleHP3NPInelasticFS, "F20/");
    clist->Register(new G4ParticleHPN2PInelasticFS, "F21/");
    clist->Register(new G4ParticleHPNPAInelasticFS, "F22/");
    clist->Register(new G4ParticleHPPInelasticFS, "F23/");
    clist->Register(new G4ParticleHPDInelasticFS, "F24/");
    clist->Register(new G4ParticleHPTInelasticFS, "F25/");
    clist->Register(new G4ParticleHPHe3InelasticFS, "F26/");
    clist->Register(new G4ParticleHPAInelasticFS, "F27/");
    clist->Register(new G4ParticleHP2AInelasticFS, "F28/");
    clist->Register(new G4ParticleHP3AInelasticFS, "F29/");
    clist->Register(new G4ParticleHP2PInelasticFS, "F30/");
    clist->Register(new G4ParticleHPPAInelasticFS, "F31/");
    clist->Register(new G4ParticleHPD2AInelasticFS, "F32/");
    clist->Register(new G4ParticleHPT2AInelasticFS, "F33/");
    clist->Register(new G4ParticleHPPDInelasticFS, "F34/");
    clist->Register(new G4ParticleHPPTInelasticFS, "F35/");
    clist->Register(new G4ParticleHPDAInelasticFS, "F36/");
#ifdef G4VERBOSE
    if (fManager->GetVerboseLevel() > 1) {
      G4cout << "ParticleHP::Inelastic for " 
	     << theProjectile->GetParticleName() << " off " 
	     << (*table)[i]->GetName() << G4endl;
    }
#endif
  }
  l.unlock();
}

void G4ParticleHPInelastic::ModelDescription(std::ostream& outFile) const
{
  outFile << "High Precision (HP) model for inelastic reaction of "
	  << theProjectile->GetParticleName() << " below 20MeV\n";
}
