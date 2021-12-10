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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file lArTPCSD.cc
/// \brief Implementation of the CaTS::lArTPCSD class

// Geant4 headers 
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#ifdef WITH_G4OPTICKS
#  include "G4Opticks.hh"
#  include "TrackInfo.hh"
#  include "OpticksGenstep.h"
#  include "OpticksFlags.hh"
#  include "G4OpticksHit.hh"
#  include "G4Cerenkov.hh"
#  include "G4Event.hh"
#  include "G4MaterialPropertiesTable.hh"
#  include "G4PhysicalConstants.hh"
#  include "G4RunManager.hh"
#  include "G4SteppingManager.hh"
#  include "G4SystemOfUnits.hh"
#  include "G4UnitsTable.hh"
#  include "G4VProcess.hh"
#  include "G4VRestDiscreteProcess.hh"
#  include "PhotonSD.hh"
#  include "G4Cerenkov.hh"
#  include "G4Scintillation.hh"
#  include "G4Version.hh"
#endif
// project headers
#include "lArTPCSD.hh"
#include "ConfigurationManager.hh"

lArTPCSD::lArTPCSD(G4String name)
  : G4VSensitiveDetector(name)
{
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  if(verbose)
  {
    G4cout << collectionName.size() << "   lArTPCSD name:  " << name
           << " collection Name: " << HCname << G4endl;
  }
  fHCID = -1;
  first = true;
}

void lArTPCSD::Initialize(G4HCofThisEvent* hce)
{
  flArTPCHitsCollection =
    new lArTPCHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
  {
    if(verbose)
    {
      G4cout << "lArTPCSD::Initialize:  " << SensitiveDetectorName << "   "
             << collectionName[0] << G4endl;
    }
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hce->AddHitsCollection(fHCID, flArTPCHitsCollection);
}

G4bool lArTPCSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return false;
  // only deal with charged particles
  G4Track* aTrack = aStep->GetTrack();
  G4double charge = aTrack->GetDynamicParticle()->GetCharge();
  if(charge == 0)
    return false;
  G4double ds       = aStep->GetStepLength();
  lArTPCHit* newHit = new lArTPCHit(
    NumElectrons(edep, ds), aStep->GetPostStepPoint()->GetPosition().getX(),
    aStep->GetPostStepPoint()->GetPosition().getY(),
    aStep->GetPostStepPoint()->GetPosition().getZ());
  flArTPCHitsCollection->insert(newHit);
#ifdef WITH_G4OPTICKS
  if(ConfigurationManager::getInstance()->isEnable_opticks())
  {
    if(first)
    {
      aMaterial     = aTrack->GetMaterial();
      materialIndex = aMaterial->GetIndex();
      if(verbose)
      {
        G4cout << "*******************************" << G4endl;
        G4cout << "RadiatorSD::ProcessHits initializing Material:  "
               << aMaterial->GetName() << " " << G4endl;
        G4cout << "RadiatorSD::ProcessHits: Name "
               << aStep->GetPreStepPoint()
                    ->GetPhysicalVolume()
                    ->GetLogicalVolume()
                    ->GetName()
               << G4endl;
      }
      aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
      if(verbose)
      {
        aMaterialPropertiesTable->DumpTable();
      }
      //
      // properties related to Scintillation
      //
#  if(G4VERSION_NUMBER > 1072)
      YieldRatio =
        aMaterialPropertiesTable->GetConstProperty(kSCINTILLATIONYIELD1) /
        aMaterialPropertiesTable->GetConstProperty(
          kSCINTILLATIONYIELD2);  // slowerRatio,
      FastTimeConstant = aMaterialPropertiesTable->GetConstProperty(
        kSCINTILLATIONTIMECONSTANT1);  // TimeConstant,
      SlowTimeConstant = aMaterialPropertiesTable->GetConstProperty(
        kSCINTILLATIONTIMECONSTANT2);  // slowerTimeConstant,
#  else
      Fast_Intensity = aMaterialPropertiesTable->GetProperty(kFASTCOMPONENT);
      Slow_Intensity = aMaterialPropertiesTable->GetProperty(kSLOWCOMPONENT);
      YieldRatio     = aMaterialPropertiesTable->GetConstProperty(kYIELDRATIO);
#  endif
      ScintillationType = Slow;
      //
      // properties related to Cerenkov
      //
      Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
#  if(G4VERSION_NUMBER > 1072)
      Pmin = Rindex->GetMinEnergy();
      Pmax = Rindex->GetMaxEnergy();
#  else
      Pmin           = Rindex->GetMinLowEdgeEnergy();
      Pmax           = Rindex->GetMaxLowEdgeEnergy();
#  endif
      dp = Pmax - Pmin;
      if(verbose)
      {
        G4cout << "nMax: " << nMax << "Pmin: " << Pmin << "Pmax: " << Pmax
               << "dp: " << dp << G4endl;
        Rindex->DumpValues();
      }
      //
      first = false;
    }
    G4int Sphotons = 0;  // number of scintillation photons this step
    G4int Cphotons = 0;  // number of Cerenkov photons this step
    //
    // info needed for generating Cerenkov photons on the GPU;
    //
    G4double maxCos                      = 0.0;
    G4double maxSin2                     = 0.0;
    G4double beta                        = 0.0;
    G4double beta1                       = 0.0;
    G4double beta2                       = 0.0;
    G4double BetaInverse                 = 0.0;
    G4double MeanNumberOfPhotons1        = 0.0;
    G4double MeanNumberOfPhotons2        = 0.0;
    G4SteppingManager* fpSteppingManager = G4EventManager::GetEventManager()
                                             ->GetTrackingManager()
                                             ->GetSteppingManager();
    G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();
    if(stepStatus != fAtRestDoItProc)
    {
      G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
      size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
      for(size_t i3 = 0; i3 < MAXofPostStepLoops; i3++)
      {
        if((*procPost)[i3]->GetProcessName() == "Cerenkov")
        {
          G4Cerenkov* proc = (G4Cerenkov*) (*procPost)[i3];
          thePhysicsTable  = proc->GetPhysicsTable();
          CerenkovAngleIntegrals =
            (G4PhysicsOrderedFreeVector*) ((*thePhysicsTable)(materialIndex));
          Cphotons = proc->GetNumPhotons();
          if(Cphotons > 0)
          {
            beta1       = aStep->GetPreStepPoint()->GetBeta();
            beta2       = aStep->GetPostStepPoint()->GetBeta();
            beta        = (beta1 + beta2) * 0.5;
            BetaInverse = 1. / beta;
            maxCos      = BetaInverse / nMax;
            maxSin2     = (1.0 - maxCos) * (1.0 + maxCos);
            MeanNumberOfPhotons1 =
              proc->GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
            MeanNumberOfPhotons2 =
              proc->GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);
          }
        }
        if((*procPost)[i3]->GetProcessName() == "Scintillation")
        {
          G4Scintillation* proc1 = (G4Scintillation*) (*procPost)[i3];
          Sphotons               = proc1->GetNumPhotons();
        }
      }
    }
    tSphotons += Sphotons;
    tCphotons += Cphotons;
    G4ThreeVector deltaPosition = aStep->GetDeltaPosition();
    G4double ScintillationTime  = 0. * ns;
    G4int scntId                = 1;
    G4StepPoint* pPreStepPoint  = aStep->GetPreStepPoint();
    G4ThreeVector x0            = pPreStepPoint->GetPosition();
    G4ThreeVector p0            = aStep->GetDeltaPosition().unit();
    //
    // harvest the Scintillation photon gensteps:
    //
    if(Sphotons > 0)
    {
      G4double ScintillationRiseTime = 0.0;
      G4Opticks::Get()->collectGenstep_G4Scintillation_1042(
        aTrack, aStep, Sphotons, scntId, ScintillationTime,
        ScintillationRiseTime);
    }
    //
    // harvest the Cerenkov photon gensteps:
    //
    if(Cphotons > 0)
    {
      G4Opticks::Get()->collectGenstep_G4Cerenkov_1042(
        aTrack, aStep, Cphotons, BetaInverse, Pmin, Pmax, maxCos, maxSin2,
        MeanNumberOfPhotons1, MeanNumberOfPhotons2);
    }
    G4Opticks* g4ok      = G4Opticks::Get();
    G4RunManager* rm     = G4RunManager::GetRunManager();
    const G4Event* event = rm->GetCurrentEvent();
    G4int eventid        = event->GetEventID();
    G4OpticksHit hit;
    unsigned num_photons = g4ok->getNumPhotons();
    if(num_photons > ConfigurationManager::getInstance()->getMaxPhotons())
    {
      g4ok->propagateOpticalPhotons(eventid);
      G4HCtable* hctable = G4SDManager::GetSDMpointer()->GetHCtable();
      for(G4int i = 0; i < hctable->entries(); ++i)
      {
        std::string sdn   = hctable->GetSDname(i);
        std::size_t found = sdn.find("Photondetector");
        if(found != std::string::npos)
        {
          PhotonSD* aSD =
            (PhotonSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector(
              sdn);
          aSD->AddOpticksHits();
        }
      }
      g4ok->reset();
    }
  }
#endif
  return true;
}

void lArTPCSD::EndOfEvent(G4HCofThisEvent*)
{
  tSphotons    = 0;
  tCphotons    = 0;
  G4int NbHits = flArTPCHitsCollection->entries();
  if(verbose)
  {
    G4cout << " Number of lArTPCHits:  " << NbHits << G4endl;
  }
}

G4double lArTPCSD::NumElectrons(G4double edep, G4double ds)
{
  // Nucl.Instrum.Meth.A523:275-286,2004
  G4double fGeVToElectrons = 4.237e+07;
  G4double fModBoxA        = 0.930;
  G4double fModBoxB        = 0.212;
  G4double EFieldStep      = 0.5;
  G4double recomb          = 0.0;
  G4double dEdx            = (ds <= 0.0) ? 0.0 : edep / ds;
  if(dEdx < 1.)
    dEdx = 1.;
  if(ds > 0)
  {
    G4double Xi = fModBoxB * dEdx / EFieldStep;
    recomb      = std::log(fModBoxA + Xi) / Xi;
  }
  else
  {
    recomb = 0.0;
  }
  G4double fNumIonElectrons = fGeVToElectrons * 1.e-3 * edep * recomb;
  return fNumIonElectrons;
}
