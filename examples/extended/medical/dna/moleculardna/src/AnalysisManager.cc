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
//
/// file:
/// brief:

#include "AnalysisManager.hh"
#include "AnalysisMessenger.hh"
#include "ChromosomeHit.hh"
#include "DNAGeometry.hh"
#include "DNAHit.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "UtilityFunctions.hh"
#include <fstream>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::AnalysisManager()
{
  fAnalysisManager    = G4AnalysisManager::Instance();
  fpAnalysisMessenger = new AnalysisMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::~AnalysisManager() { delete fpAnalysisMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Initialize()
{
  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
  if(run->GetRunID() == 0)
  {
    G4cout << "AnalysisManager::Initialize() GetRunID : " << run->GetRunID()
           << G4endl;
    fAnalysisManager->SetDefaultFileType("root");
    fAnalysisManager->SetVerboseLevel(1);

    fAnalysisManager->SetNtupleDirectoryName("tuples");
    fAnalysisManager->SetHistoDirectoryName("hists");

    fAnalysisManager->CreateH1("ssb_counts", "SSBs", 1000, 0, 1000);
    fAnalysisManager->CreateH1("ssb_energies_ev", "SSB Energies", 500, 0, 500);
    // fAnalysisManager->CreateH1("fragments", "Fragment sizes", 500, 0,
    // 500);//ORG
    fAnalysisManager->CreateH1("fragments", "Fragment sizes", 10000, 0,
                               10000);  // dousatsu
    fAnalysisManager->CreateH3("local_pos", "Strand Interaction Positions", 50,
                               -25, 25, 50, -25, 25, 50, -25, 25);
    fAnalysisManager->CreateH3("global_pos", "Strand Interaction Positions", 50,
                               -100, 100, 50, -100, 100, 50, -100, 100);
    fAnalysisManager->CreateH1("e_cell_kev", "Energy Deposits in Cell", 2500, 0,
                               2500);  // dousatsu

    fAnalysisManager->CreateNtuple("primary_source", "Primary Source");
    fAnalysisManager->CreateNtupleSColumn("Primary");
    fAnalysisManager->CreateNtupleDColumn("Energy");
    fAnalysisManager->CreateNtupleDColumn("PosX_um");
    fAnalysisManager->CreateNtupleDColumn("PosY_um");
    fAnalysisManager->CreateNtupleDColumn("PosZ_um");
    fAnalysisManager->CreateNtupleDColumn("MomX");
    fAnalysisManager->CreateNtupleDColumn("MomY");
    fAnalysisManager->CreateNtupleDColumn("MomZ");
    fAnalysisManager->CreateNtupleDColumn("StopPosX_um");
    fAnalysisManager->CreateNtupleDColumn("StopPosY_um");
    fAnalysisManager->CreateNtupleDColumn("StopPosZ_um");
    fAnalysisManager->CreateNtupleDColumn("TraLen_cell_um");
    fAnalysisManager->CreateNtupleDColumn("TraLen_chro_um");
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("chromosome_hits",
                                   "Energy Deposits in Chromosomes");
    fAnalysisManager->CreateNtupleSColumn("chromosome");
    fAnalysisManager->CreateNtupleDColumn("e_chromosome_kev");
    fAnalysisManager->CreateNtupleDColumn("e_dna_kev");
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("classification",
                                   "Break Complexity Frequency");
    fAnalysisManager->CreateNtupleSColumn("Primary");
    fAnalysisManager->CreateNtupleDColumn("Energy");
    fAnalysisManager->CreateNtupleIColumn("None");
    fAnalysisManager->CreateNtupleIColumn("SSB");
    fAnalysisManager->CreateNtupleIColumn("SSBp");
    fAnalysisManager->CreateNtupleIColumn("2SSB");
    fAnalysisManager->CreateNtupleIColumn("DSB");
    fAnalysisManager->CreateNtupleIColumn("DSBp");
    fAnalysisManager->CreateNtupleIColumn("DSBpp");
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("source", "Break Source Frequency");
    fAnalysisManager->CreateNtupleSColumn("Primary");
    fAnalysisManager->CreateNtupleDColumn("Energy");
    fAnalysisManager->CreateNtupleIColumn("None");
    fAnalysisManager->CreateNtupleIColumn("SSBd");
    fAnalysisManager->CreateNtupleIColumn("SSBi");
    fAnalysisManager->CreateNtupleIColumn("SSBm");
    fAnalysisManager->CreateNtupleIColumn("DSBd");
    fAnalysisManager->CreateNtupleIColumn("DSBi");
    fAnalysisManager->CreateNtupleIColumn("DSBm");
    fAnalysisManager->CreateNtupleIColumn("DSBh");
    fAnalysisManager->FinishNtuple();

    fAnalysisManager->CreateNtuple("damage", "DNA damage locations");
    fAnalysisManager->CreateNtupleIColumn("Event");
    fAnalysisManager->CreateNtupleSColumn("Primary");
    fAnalysisManager->CreateNtupleDColumn("Energy");
    fAnalysisManager->CreateNtupleSColumn("TypeClassification");
    fAnalysisManager->CreateNtupleSColumn("SourceClassification");
    fAnalysisManager->CreateNtupleDColumn("Position_x_um");
    fAnalysisManager->CreateNtupleDColumn("Position_y_um");
    fAnalysisManager->CreateNtupleDColumn("Position_z_um");
    fAnalysisManager->CreateNtupleDColumn("Size_nm");
    fAnalysisManager->CreateNtupleIColumn("FragmentLength");
    fAnalysisManager->CreateNtupleIColumn("BaseDamage");
    fAnalysisManager->CreateNtupleIColumn("StrandDamage");
    fAnalysisManager->CreateNtupleIColumn("DirectBreaks");
    fAnalysisManager->CreateNtupleIColumn("IndirectBreaks");
    fAnalysisManager->CreateNtupleIColumn("EaqBaseHits");
    fAnalysisManager->CreateNtupleIColumn("EaqStrandHits");
    fAnalysisManager->CreateNtupleIColumn("OHBaseHits");
    fAnalysisManager->CreateNtupleIColumn("OHStrandHits");
    fAnalysisManager->CreateNtupleIColumn("HBaseHits");
    fAnalysisManager->CreateNtupleIColumn("HStrandHits");
    fAnalysisManager->CreateNtupleDColumn("EnergyDeposited_eV");
    fAnalysisManager->CreateNtupleIColumn("InducedBreaks");
    fAnalysisManager->CreateNtupleIColumn("Chain");
    fAnalysisManager->CreateNtupleIColumn("Strand");    // WG
    fAnalysisManager->CreateNtupleDColumn("BasePair");  // dousatsu
    fAnalysisManager->CreateNtupleSColumn("Name");
    fAnalysisManager->FinishNtuple();
  }

  if(!fFileName)
  {
    fFileName = "molecular-dna";
  }
  fAnalysisManager->OpenFile(fFileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Close()
{
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::ProcessDNAHitsVector(
  const std::vector<const DNAHit*>& dnaHits)
{
  // Check we have DNA Geometry (runs once)
  if(fpDNAGeometry == nullptr)
  {
    auto det = dynamic_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fpDNAGeometry = det->GetDNAGeometry();
  }

  if(dnaHits.empty())
  {
    return;
  }
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();
  // SSBs
  for(auto dnaHit : dnaHits)
  {
    fAnalysisManager->FillH1(fAnalysisManager->GetH1Id("ssb_energies_ev"),
                             dnaHit->GetEnergy() / eV);
    if(fChainToSave == -1)
    {
      // save all chains
      fAnalysisManager->FillH3(fAnalysisManager->GetH3Id("local_pos"),
                               dnaHit->GetLocalPosition().getX() / nm,
                               dnaHit->GetLocalPosition().getY() / nm,
                               dnaHit->GetLocalPosition().getZ() / nm);
      fAnalysisManager->FillH3(fAnalysisManager->GetH3Id("global_pos"),
                               dnaHit->GetPosition().getX() / nm,
                               dnaHit->GetPosition().getY() / nm,
                               dnaHit->GetPosition().getZ() / nm);
    }
    else if(fChainToSave == dnaHit->GetChainIdx())
    {  // what is different ?
      // save only specified chain
      fAnalysisManager->FillH3(fAnalysisManager->GetH3Id("local_pos"),
                               dnaHit->GetLocalPosition().getX() / nm,
                               dnaHit->GetLocalPosition().getY() / nm,
                               dnaHit->GetLocalPosition().getZ() / nm);
      fAnalysisManager->FillH3(fAnalysisManager->GetH3Id("global_pos"),
                               dnaHit->GetPosition().getX() / nm,
                               dnaHit->GetPosition().getY() / nm,
                               dnaHit->GetPosition().getZ() / nm);
    }
  }
  fAnalysisManager->FillH1(fAnalysisManager->GetH1Id("ssb_counts"),
                           dnaHits.size());

  // DSBs
  // Sorting
  // Make a map relating chromosome/chain to its binary tree
  std::map<std::pair<G4String, G4int>, BinaryTree*> treemap;

  // Populate the binary tree
  // G4cout << "Building Binary Trees for " << dnaHits.size() << " Nodes" <<
  // G4endl;
  for(auto it = dnaHits.begin(); it != dnaHits.end(); ++it)
  {
    std::pair<G4String, G4int> key =
      std::make_pair((*it)->GetChromosome(), (*it)->GetChainIdx());
    if(treemap.find(key) == treemap.end())
    {
      treemap[key] = new BinaryTree();
      if(fAnalysisManager->GetVerboseLevel() >= 2)
      {
        G4cout << "Constructing binary tree. Chromosome: " << key.first
               << " Chain: " << key.second << G4endl;
      }
    }
    treemap.at(key)->Insert(*it);
    if(fAnalysisManager->GetVerboseLevel() >= 2)
    {
      G4cout << "Adding hit to binary tree. Chromosome: " << key.first
             << " Chain: " << key.second << G4endl;
    }
  }

  // Analysis
  // Create a vector to stock damage records for the event
  // G4cout << "Building Damage Records" << G4endl;
  std::vector<DamageRecord*> damageRecords;
  for(auto& it : treemap)
  {
    if(fAnalysisManager->GetVerboseLevel() >= 2)
    {
      G4cout << "Analysing hits for chromosome: " << it.first.first
             << " Chain: " << it.first.second << G4endl;
    }
    DNAHit* currentHit = it.second->First();
    DNAHit* nextHit    = currentHit;

    // Runs while their are still hits
    while(nextHit)
    {
      // Create record for this fragment
      if(currentHit->GetBasePairIdx() > 9223372036854775807)
      {
        G4cout << " SEE AnalysisManager !!!" << G4endl;
        abort();
      }  // dousatsu
      // if (currentHit->GetBasePairIdx()>2147483647){G4cout<<" SEE
      // AnalysisManager !!!"<<G4endl;abort();}//ORG
      auto idx         = (int64_t) currentHit->GetBasePairIdx();  // dousatsu
      int64_t startidx = (idx > 10) ? (idx - 10) : 0;             // dousatsu
      // G4int idx = currentHit->GetBasePairIdx();//ORG
      // G4int startidx = (idx > 10) ? (idx - 10) : 0;//ORG
      G4String name = it.first.first + "_" + std::to_string(it.first.second) +
                      "_" + std::to_string(startidx);  // check!!!!
      auto* record = new DamageRecord(name, idx, currentHit->GetPlacementIdx(),
                                      currentHit->GetChainIdx());
      damageRecords.push_back(record);
      BasePairDamageRecord* bp;
      int64_t gap = 0;  // dousatsu

      // bookend with 10 bps on either side of no damage
      record->AddEmptyBPDamage(idx - startidx);
      // continues until fragment ends
      while(true)
      {
        bp                      = new BasePairDamageRecord;
        bp->fStrand1Energy      = currentHit->GetStrand1Energy();
        bp->fStrand2Energy      = currentHit->GetStrand2Energy();
        bp->fBp1Energy          = currentHit->GetBP1Energy();
        bp->fBp2Energy          = currentHit->GetBP2Energy();
        bp->fBp1IndirectEvt     = (currentHit->GetBase1Rad() != nullptr);
        bp->fBp2IndirectEvt     = (currentHit->GetBase2Rad() != nullptr);
        bp->fStrand1IndirectEvt = (currentHit->GetStrand1Rad() != nullptr);
        bp->fStrand2IndirectEvt = (currentHit->GetStrand2Rad() != nullptr);
        bp->fbp1DirectDmg       = false;  // No physical damage on BPs
        bp->fbp2DirectDmg       = false;  // No physical damage on BPs

        bp->fStrand1DirectDmg =
          fpDNAGeometry->GetDamageModel()->IsDirectStrandBreak(
            bp->fStrand1Energy);
        bp->fStrand2DirectDmg =
          fpDNAGeometry->GetDamageModel()->IsDirectStrandBreak(
            bp->fStrand2Energy);

        //        if(bp->fStrand1DirectDmg) bp->fStrandIdx = 1;
        //        else if(bp->fStrand2DirectDmg) strandID = 2;

        if(bp->fBp1IndirectEvt)
        {
          bp->fBp1IndirectDmg =
            fpDNAGeometry->GetDamageModel()->IsIndirectBaseDamage(
              currentHit->GetBase1Rad());
          if(bp->fBp1IndirectDmg)
          {
            bp->fbp1InducedBreak =
              fpDNAGeometry->GetDamageModel()->IsInducedStrandBreak(
                currentHit->GetBase1Rad());
          }
        }

        if(bp->fBp2IndirectEvt)
        {
          bp->fBp2IndirectDmg =
            fpDNAGeometry->GetDamageModel()->IsIndirectBaseDamage(
              currentHit->GetBase2Rad());
          if(bp->fBp2IndirectDmg)
          {
            bp->fbp2InducedBreak =
              fpDNAGeometry->GetDamageModel()->IsInducedStrandBreak(
                currentHit->GetBase2Rad());
          }
        }

        if(bp->fStrand1IndirectEvt)
        {
          bp->fStrand1IndirectDmg =
            fpDNAGeometry->GetDamageModel()->IsIndirectStrandDamage(
              currentHit->GetStrand1Rad());
        }
        if(bp->fStrand2IndirectEvt)
        {
          bp->fStrand2IndirectDmg =
            fpDNAGeometry->GetDamageModel()->IsIndirectStrandDamage(
              currentHit->GetStrand2Rad());
        }
        // Record radical-base/strand damage.
        // it is possible for base/strand damage to not correspond to
        // direct hits due to the damage induction probabilities.
        // This mainly affects strands.
        if(bp->fBp1IndirectDmg)
        {
          record->AddBaseHit(currentHit->GetBase1Rad()->GetDefinition());
        }
        if(bp->fBp2IndirectDmg)
        {
          record->AddBaseHit(currentHit->GetBase2Rad()->GetDefinition());
        }
        if(bp->fStrand1IndirectDmg)
        {
          record->AddStrandHit(currentHit->GetStrand1Rad()->GetDefinition());
          //         strandID = 1;
        }
        if(bp->fStrand2IndirectDmg)
        {
          record->AddStrandHit(currentHit->GetStrand2Rad()->GetDefinition());
          //          strandID = 2;
        }

        record->AddBasePairDamage(bp, currentHit->GetPosition());

        nextHit = it.second->Next(currentHit);
        if(nextHit == nullptr)
        {
          break;
        }

        gap = nextHit->GetBasePairIdx() - currentHit->GetBasePairIdx();
        if(fFragmentGap > 0)
        {
          // case 1: continuous strand
          if(gap > fFragmentGap)
          {
            currentHit = nextHit;
            break;
          }
          else
          {
            record->AddEmptyBPDamage(gap);
          }
        }
        else if(fFragmentGap == 0)
        {
          // case 2: individual placements, not joined
          // ie. plasmids etc. Each placement is a separate strand
          if(currentHit->GetPlacementIdx() != nextHit->GetPlacementIdx())
          {
            if(fpDNAGeometry->GetVerbosity() > 1)
            {
              G4cout << "Analysis passing to a new placement" << G4endl;
            }
            currentHit = nextHit;
            break;
          }
          else
          {
            record->AddEmptyBPDamage(gap);
          }
        }
        else
        {
          G4Exception("MolecularAnalaysisManager", "ERR_FRAGMENT_VALUE",
                      FatalException,
                      "The value set in /analysisDNA/fragmentGap is bad");
        }
        currentHit = nextHit;
      }
      record->AddEmptyBPDamage(10);
      // end bookend
    }
  }
  // Count the number of breaks and print out the records
  // std::map<G4String, G4int> breakmap;
  std::map<complexityEnum, G4int> breakmap;

  breakmap[DSBplusplus]    = 0;
  breakmap[DSBplus]        = 0;
  breakmap[DSB]            = 0;
  breakmap[twoSSB]         = 0;
  breakmap[SSBplus]        = 0;
  breakmap[SSB]            = 0;
  breakmap[NoneComplexity] = 0;
  std::map<sourceEnum, G4int> sourcemap;
  sourcemap[SSBd]      = 0;
  sourcemap[SSBi]      = 0;
  sourcemap[SSBm]      = 0;
  sourcemap[DSBh]      = 0;
  sourcemap[DSBm]      = 0;
  sourcemap[DSBd]      = 0;
  sourcemap[DSBi]      = 0;
  sourcemap[undefined] = 0;

  for(auto& damageRecord : damageRecords)
  {
    DamageClassification* classif =
      damageRecord->GetClassification(fDSBDistance);
    breakmap[classif->fComplexity]++;
    sourcemap[classif->fSource]++;
    fAnalysisManager->FillH1(fAnalysisManager->GetH1Id("fragments"),
                             damageRecord->GetSize());

    ///////////dousatsu--------------------------------------------
    fAnalysisManager->FillNtupleIColumn(4, 0, event->GetEventID());
    fAnalysisManager->FillNtupleSColumn(4, 1,
                                        event->GetPrimaryVertex()
                                          ->GetPrimary()
                                          ->GetParticleDefinition()
                                          ->GetParticleName());
    fAnalysisManager->FillNtupleDColumn(
      4, 2, event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());

    G4String complexityString = "None";
    switch(classif->fComplexity)
    {
      case SSB:
        complexityString = "SSB";
        break;
      case SSBplus:
        complexityString = "SSB+";
        break;
      case twoSSB:
        complexityString = "2SSB";
        break;
      case DSB:
        complexityString = "DSB";
        break;
      case DSBplus:
        complexityString = "DSB+";
        break;
      case DSBplusplus:
        complexityString = "DSB++";
        break;
      default:
        complexityString = "None";
        break;
    }
    fAnalysisManager->FillNtupleSColumn(4, 3, complexityString);

    G4String sourceString = "undefined";
    switch(classif->fSource)
    {
      case SSBd:
        sourceString = "SSBd";
        break;
      case SSBi:
        sourceString = "SSBi";
        break;
      case SSBm:
        sourceString = "SSBm";
        break;
      case DSBh:
        sourceString = "DSBh";
        break;
      case DSBm:
        sourceString = "DSBm";
        break;
      case DSBd:
        sourceString = "DSBd";
        break;
      case DSBi:
        sourceString = "DSBi";
        break;
      default:
        sourceString = "undefined";
        break;
    }

    fAnalysisManager->FillNtupleSColumn(4, 4, sourceString);
    fAnalysisManager->FillNtupleDColumn(
      4, 5, damageRecord->GetMeanPosition().getX() / um);
    fAnalysisManager->FillNtupleDColumn(
      4, 6, damageRecord->GetMeanPosition().getY() / um);
    fAnalysisManager->FillNtupleDColumn(
      4, 7, damageRecord->GetMeanPosition().getZ() / um);
    fAnalysisManager->FillNtupleDColumn(4, 8,
                                        damageRecord->GetMeanDistance() / nm);
    fAnalysisManager->FillNtupleIColumn(4, 9, damageRecord->GetSize());
    fAnalysisManager->FillNtupleIColumn(4, 10, classif->fbaseDmg);
    fAnalysisManager->FillNtupleIColumn(4, 11, classif->fStrandDmg);
    fAnalysisManager->FillNtupleIColumn(4, 12, classif->fDirectBreaks);
    fAnalysisManager->FillNtupleIColumn(4, 13, classif->fIndirectBreaks);
    fAnalysisManager->FillNtupleIColumn(4, 14, damageRecord->GetEaqBaseHits());
    fAnalysisManager->FillNtupleIColumn(4, 15,
                                        damageRecord->GetEaqStrandHits());
    fAnalysisManager->FillNtupleIColumn(4, 16, damageRecord->GetOHBaseHits());
    fAnalysisManager->FillNtupleIColumn(4, 17, damageRecord->GetOHStrandHits());
    fAnalysisManager->FillNtupleIColumn(4, 18, damageRecord->GetHBaseHits());
    fAnalysisManager->FillNtupleIColumn(4, 19, damageRecord->GetHStrandHits());
    fAnalysisManager->FillNtupleDColumn(4, 20, damageRecord->GetEnergy() / eV);
    fAnalysisManager->FillNtupleIColumn(4, 21, classif->fInducedBreaks);
    fAnalysisManager->FillNtupleIColumn(4, 22, damageRecord->GetChainIdx());
    fAnalysisManager->FillNtupleIColumn(4, 23, damageRecord->GetPlacementIdx());
    fAnalysisManager->FillNtupleDColumn(
      4, 24, (G4double) damageRecord->GetStartBPIdx());
    // fAnalysisManager->FillNtupleIColumn(4, 23, (*it)->GetStartBPIdx());//ORG
    fAnalysisManager->FillNtupleSColumn(4, 25, damageRecord->GetName());
    fAnalysisManager->AddNtupleRow(4);
    delete classif;
    //
    // TODO: Use (*it)->GetPlacementIdx() to work out if the damage
    // is continuously joined to a past event or not.
    // If yes, save the fragment length.
    //
    if(fSaveStrands)
    {
      damageRecord->PrintRecord(fStrandDirectory + "/" +
                                std::to_string(event->GetEventID()) + "_" +
                                damageRecord->GetName() + ".txt");
    }
    if(fAnalysisManager->GetVerboseLevel() >= 2)
      damageRecord->PrintRecord("");
  }

  fAnalysisManager->FillNtupleSColumn(2, 0,
                                      event->GetPrimaryVertex()
                                        ->GetPrimary()
                                        ->GetParticleDefinition()
                                        ->GetParticleName());
  fAnalysisManager->FillNtupleDColumn(
    2, 1, event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());
  fAnalysisManager->FillNtupleIColumn(2, 2, breakmap[NoneComplexity]);
  fAnalysisManager->FillNtupleIColumn(2, 3, breakmap[SSB]);
  fAnalysisManager->FillNtupleIColumn(2, 4, breakmap[SSBplus]);
  fAnalysisManager->FillNtupleIColumn(2, 5, breakmap[twoSSB]);
  fAnalysisManager->FillNtupleIColumn(2, 6, breakmap[DSB]);
  fAnalysisManager->FillNtupleIColumn(2, 7, breakmap[DSBplus]);
  fAnalysisManager->FillNtupleIColumn(2, 8, breakmap[DSBplusplus]);
  fAnalysisManager->AddNtupleRow(2);

  fAnalysisManager->FillNtupleSColumn(3, 0,
                                      event->GetPrimaryVertex()
                                        ->GetPrimary()
                                        ->GetParticleDefinition()
                                        ->GetParticleName());
  fAnalysisManager->FillNtupleDColumn(
    3, 1, event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());
  fAnalysisManager->FillNtupleIColumn(3, 2, sourcemap[undefined]);
  fAnalysisManager->FillNtupleIColumn(3, 3, sourcemap[SSBd]);
  fAnalysisManager->FillNtupleIColumn(3, 4, sourcemap[SSBi]);
  fAnalysisManager->FillNtupleIColumn(3, 5, sourcemap[SSBm]);
  fAnalysisManager->FillNtupleIColumn(3, 6, sourcemap[DSBd]);
  fAnalysisManager->FillNtupleIColumn(3, 7, sourcemap[DSBi]);
  fAnalysisManager->FillNtupleIColumn(3, 8, sourcemap[DSBm]);
  fAnalysisManager->FillNtupleIColumn(3, 9, sourcemap[DSBh]);
  fAnalysisManager->AddNtupleRow(3);

  for(auto it : damageRecords)
  {
    delete it;
  }

  // Cleanup, delete binary trees
  for(const auto& it : treemap)
  {
    delete it.second;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::ProcessChromosomeHitMap(
  const std::map<uint32_t, ChromosomeHit*>& chromosomes)
{
  for(auto chromosome : chromosomes)
  {
    fAnalysisManager->FillNtupleSColumn(1, 0, chromosome.second->GetName());
    fAnalysisManager->FillNtupleDColumn(
      1, 1, chromosome.second->GetChromosomeEdep() / keV);
    fAnalysisManager->FillNtupleDColumn(1, 2,
                                        chromosome.second->GetDNAEdep() / keV);
    fAnalysisManager->AddNtupleRow(1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::ProcessPrimary(const G4ThreeVector& primstoppos,
                                     const G4double& tlcell,
                                     const G4double& tlchro)
{
  const G4Event* event =
    G4EventManager::GetEventManager()->GetConstCurrentEvent();
  fAnalysisManager->FillNtupleSColumn(0, 0,
                                      event->GetPrimaryVertex()
                                        ->GetPrimary()
                                        ->GetParticleDefinition()
                                        ->GetParticleName());
  fAnalysisManager->FillNtupleDColumn(
    0, 1, event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());
  fAnalysisManager->FillNtupleDColumn(0, 2,
                                      event->GetPrimaryVertex()->GetX0() / um);
  fAnalysisManager->FillNtupleDColumn(0, 3,
                                      event->GetPrimaryVertex()->GetY0() / um);
  fAnalysisManager->FillNtupleDColumn(0, 4,
                                      event->GetPrimaryVertex()->GetZ0() / um);

  G4ThreeVector mom =
    event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection();
  fAnalysisManager->FillNtupleDColumn(0, 5, mom.x() / mom.mag());
  fAnalysisManager->FillNtupleDColumn(0, 6, mom.y() / mom.mag());
  fAnalysisManager->FillNtupleDColumn(0, 7, mom.z() / mom.mag());
  fAnalysisManager->FillNtupleDColumn(0, 8, primstoppos.x() / um);
  fAnalysisManager->FillNtupleDColumn(0, 9, primstoppos.y() / um);
  fAnalysisManager->FillNtupleDColumn(0, 10, primstoppos.z() / um);
  fAnalysisManager->FillNtupleDColumn(0, 11, tlcell / um);
  fAnalysisManager->FillNtupleDColumn(0, 12, tlchro / um);
  fAnalysisManager->AddNtupleRow(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::ProcessCellEdep(const G4double& edep)
{
  fAnalysisManager->FillH1(fAnalysisManager->GetH1Id("e_cell_kev"), edep / keV);
}

////////////////////////////////////////////////////////////////////////////////
////// Binary Tree Methods
////////////////////////////////////////////////////////////////////////////////

BinaryTree::BinaryTree() { fRoot = nullptr; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BinaryTree::~BinaryTree() { Destroy_tree(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Makes own internal copy of hit
void BinaryTree::Insert(const DNAHit* hit)
{
  auto newHit = new DNAHit(*hit);
  if(fRoot == nullptr)
  {
    Node* node    = new Node;
    node->fleft   = nullptr;
    node->fright  = nullptr;
    node->fparent = nullptr;
    node->fdata   = newHit;
    node->fkey    = newHit->GetBasePairIdx();
    fRoot         = node;
  }
  else
  {
    Insert_(newHit, fRoot);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit* BinaryTree::Search(const int64_t& index)  // dousatsu
{
  return Search_(index, fRoot);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BinaryTree::Destroy_tree() { Destroy_tree_(fRoot); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BinaryTree::Destroy_tree_(Node* node)
{
  if(node != nullptr)
  {
    Destroy_tree_(node->fleft);
    Destroy_tree_(node->fright);
    delete node->fdata;
    delete node;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BinaryTree::Insert_(DNAHit* hit, Node* node)
{
  if(hit->GetBasePairIdx() > 9223372036854775807)
  {
    G4cout << " SEE AnalysisManager !!!" << G4endl;
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4cout << "!!! aborting ... " << G4endl;
    G4cout << "This pair ID is so big " << hit->GetBasePairIdx() << G4endl;
    G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4Exception("BinaryTree"
                "::insert()",
                "insert000", FatalException, exceptionDescription);
  }
  int64_t key = hit->GetBasePairIdx();  // dousatsu
  // G4int key = hit->GetBasePairIdx();//ORG
  if(key < node->fkey)
  {
    if(node->fleft != nullptr)
    {
      Insert_(hit, node->fleft);
    }
    else
    {
      Node* daughter    = new Node;
      daughter->fparent = node;
      daughter->fleft   = nullptr;
      daughter->fright  = nullptr;
      daughter->fdata   = hit;
      daughter->fkey    = hit->GetBasePairIdx();
      node->fleft       = daughter;
    }
  }
  else if(key > node->fkey)
  {
    if(node->fright != nullptr)
    {
      Insert_(hit, node->fright);
    }
    else
    {
      Node* daughter    = new Node;
      daughter->fparent = node;
      daughter->fleft   = nullptr;
      daughter->fright  = nullptr;
      daughter->fdata   = hit;
      daughter->fkey    = hit->GetBasePairIdx();
      node->fright      = daughter;
    }
  }
  else
  {
    node->fdata->AddHit(*hit);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit* BinaryTree::Search_(const int64_t& key, Node* node)  // dousatsu check
{
  if(node != nullptr)
  {
    if(key < node->fkey)
    {
      return Search_(key, node->fleft);
    }
    else if(key > node->fkey)
    {
      return Search_(key, node->fright);
    }
    else
    {
      return node->fdata;
    }
  }
  else
  {
    return nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit* BinaryTree::First() const
{
  if(fRoot == nullptr)
  {
    return nullptr;
  }
  else
  {
    return First_(fRoot);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit* BinaryTree::First_(Node* node) const
{
  if(node->fleft == nullptr)
  {
    return node->fdata;
  }
  else
  {
    return First_(node->fleft);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit* BinaryTree::Next(const DNAHit* hit) const
{
  if(hit->GetBasePairIdx() > 9223372036854775807)
  {
    G4cout << " SEE AnalysisManager !!!" << G4endl;
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4cout << "!!! aborting ... " << G4endl;
    G4cout << "This pair ID is so big " << hit->GetBasePairIdx() << G4endl;
    G4cout << "!!!!!!!!!!!!!!!!!!!!" << G4endl;
    G4Exception("BinaryTree"
                "::insert()",
                "insert000", FatalException, exceptionDescription);
  }
  int64_t key = hit->GetBasePairIdx();  // dousatsu
  if(fRoot == nullptr)
  {
    return nullptr;
  }
  else
  {
    return Next_(key, fRoot);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Looks more complicated than it is
// The interface to these functions with integer keys rather than node
// objects makes this look more complicated than standard algorithms
DNAHit* BinaryTree::Next_(const int64_t& key, Node* node) const  // dousatsu
{
  if(key < node->fkey)
  {
    if(node->fleft != nullptr)
    {
      return Next_(key, node->fleft);
    }
    else  // left is NULL
    {
      return node->fdata;
    }
  }
  else  // (key >= node->key)
  {
    if(node->fright != nullptr)
    {
      return Next_(key, node->fright);
    }
    else  // right is NULL; solution is higher in tree
    {
      while(true)
      {
        node = node->fparent;
        if(node == nullptr)
        {
          return nullptr;
        }
        if(key < node->fkey)
        {
          return node->fdata;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

////////////////////////////////////////////////////////////////////////////////
////////// Damage Record
////////////////////////////////////////////////////////////////////////////////

const char* DamageRecord::fDirectDamageChar   = "D";
const char* DamageRecord::fIndirectDamageChar = "I";
const char* DamageRecord::fHitNoDamageChar    = "~";
const char* DamageRecord::fNotHitChar         = "-";
const char* DamageRecord::fBothDamageChar     = "X";

DamageRecord::DamageRecord(G4String name, const int64_t& startIndex,
                           const G4int& place_idx, const G4int& chain_idx)
  : fName(std::move(name))
  , fStartIndex(startIndex)
  , fStartPlacement(place_idx)
  , fChainIdx(chain_idx)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DamageRecord::~DamageRecord()
{
  for(auto& fDamageRecord : fDamageRecords)
  {
    delete fDamageRecord;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageRecord::PrintRecord(const G4String& filename,
                               const G4double& dsbDistance)
{
  std::stringstream strand1;
  std::stringstream bp1;
  std::stringstream bp2;
  std::stringstream strand2;
  for(auto& fDamageRecord : fDamageRecords)
  {
    strand1 << GetChar(fDamageRecord->fStrand1DirectDmg,
                       fDamageRecord->fStrand1IndirectDmg,
                       fDamageRecord->fStrand1Energy);
    bp1 << GetChar(fDamageRecord->fbp1DirectDmg, fDamageRecord->fBp1IndirectDmg,
                   fDamageRecord->fBp1Energy);
    bp2 << GetChar(fDamageRecord->fbp2DirectDmg, fDamageRecord->fBp2IndirectDmg,
                   fDamageRecord->fBp2Energy);
    strand2 << GetChar(fDamageRecord->fStrand2DirectDmg,
                       fDamageRecord->fStrand2IndirectDmg,
                       fDamageRecord->fStrand2Energy);
  }

  // print to with no filename
  if(filename.empty())
  {
    DamageClassification* classification = this->GetClassification(dsbDistance);
    G4cout << "Sequences starts at base pair " << fStartIndex
           << " in placement " << fStartPlacement << ", Chain " << fChainIdx
           << ". Total length: " << fDamageRecords.size()
           << ". Classification: " << classification->fComplexity << " "
           << classification->fSource << G4endl;
    G4cout << strand1.str() << G4endl;
    G4cout << bp1.str() << G4endl;
    G4cout << bp2.str() << G4endl;
    G4cout << strand2.str() << G4endl;
    return;
  }
  // To Text output
  std::fstream fs(filename,
                  std::fstream::in | std::fstream::out | std::fstream::trunc);
  if(fs.fail())
  {
    G4cout << "Could not open filestream" << G4endl;
  }
  DamageClassification* classification = this->GetClassification(dsbDistance);
  fs << "Sequences starts at base pair " << fStartIndex << " in placement "
     << fStartPlacement << ", Chain " << fChainIdx
     << ". Total length: " << fDamageRecords.size()
     << ". Classification: " << classification->fComplexity << " "
     << classification->fSource << std::endl;
  fs << strand1.str() << std::endl;
  fs << bp1.str() << std::endl;
  fs << bp2.str() << std::endl;
  fs << strand2.str() << std::endl;
  fs.close();
  delete classification;
}

const char* DamageRecord::GetChar(const G4bool& direct, const G4bool& indirect,
                                  const G4double& e)
{
  if(direct && indirect)
  {
    return fBothDamageChar;
  }
  else if(direct)
  {
    return fDirectDamageChar;
  }
  else if(indirect)
  {
    return fIndirectDamageChar;
  }
  else if(e > 0)
  {
    return fHitNoDamageChar;
  }
  else
  {
    return fNotHitChar;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DamageClassification* DamageRecord::GetClassification(
  const G4double& dsbDistance)
{
  /* Explanation of class because we do some sneaky bit-level manipulation to
     track breaks along the strand.

     Things to note straight off
     1) This routine makes no promise to count SSBs and DSBs correctly
      once they cease being relevant to the break complexity. That is, once
      a DSB has occurred, the SSB count becomes irrelevant
     2) The previous ten base pairs are kept in the ints lastTenStrandX
      Each bit in the number represents one base pair precedent. The bit
      is set if the strand was broken.
     3) The countX variables are used to track the number of breaks in the
      preceeding 10 base pairs on each strand. This is needed to calculate
      DSB+ events
     4) The lambda function count_bytes is an algorithm to count the number
      of set bytes in constant time for 32-bit integers. It yields for
      us the count of breaks in the previous ten bases.

     Classification of complexity:
     DSB++  2DSBs in fragment
     DSB+   1 DSB and then an SSB within 10 base pairs
     DSB  Each strand hit once within <=10 base pairs
     2SSB   Each strand hit once with  >10 base pairs separation
     SSB+   One strand hit twice or more
     SSB  One strand hit only once
     None   Nothing hit

     Classification by Source

  */
  // lambda function to count the number of bytes set in a 32-bit uint
  // https://blogs.msdn.microsoft.com/jeuge/2005/06/08/bit-fiddling-3/
  // Note: integers beginning with 0 are entered as octal!
  auto count_bytes = [](uint32_t u) {
    uint32_t uCount = u - ((u >> 1) & 033333333333) - ((u >> 2) & 011111111111);
    return ((uCount + (uCount >> 3)) & 030707070707) % 63;
  };

  auto classification = new DamageClassification();

  G4int nDSB = 0;
  G4int nSSB = 0;
  // Source classification
  G4int oldnDSB = 0;
  G4int nDSBm   = 0;
  G4int nDSBi   = 0;
  G4int nDSBd   = 0;
  G4int nDSBh   = 0;

  // Complexity Classification
  G4int nDSBPlus = 0;
  G4bool SSB1    = false;
  G4bool SSB2    = false;
  // Counters for the last ten bp
  uint32_t lastTenDirectStrand1   = 0;
  uint32_t lastTenDirectStrand2   = 0;
  uint32_t lastTenIndirectStrand1 = 0;
  uint32_t lastTenIndirectStrand2 = 0;
  uint32_t lastTenStrand1         = 0;
  uint32_t lastTenStrand2         = 0;
  uint32_t lastTenTracked1 = 0;  // We remove entries from these if they have
  uint32_t lastTenTracked2 = 0;  // been counted in DSBs
  uint32_t count1          = 0;
  uint32_t count2          = 0;
  G4bool strand1Indirect   = false;
  G4bool strand1Direct     = false;
  G4bool strand2Indirect   = false;
  G4bool strand2Direct     = false;
  G4int baseDamage         = 0;
  G4int strandDamage       = 0;
  uint32_t truncator       = std::pow(2, dsbDistance + 1) - 1;
  for(G4int ii = 0; ii != (G4int) fDamageRecords.size(); ii++)
  {
    lastTenDirectStrand1   = (lastTenDirectStrand1 << 1) & truncator;
    lastTenDirectStrand2   = (lastTenDirectStrand2 << 1) & truncator;
    lastTenIndirectStrand1 = (lastTenIndirectStrand1 << 1) & truncator;
    lastTenIndirectStrand2 = (lastTenIndirectStrand2 << 1) & truncator;

    lastTenStrand1  = lastTenStrand1 << 1;
    lastTenStrand2  = lastTenStrand2 << 1;
    lastTenTracked1 = lastTenTracked1 << 1;
    lastTenTracked2 = lastTenTracked2 << 1;
    // DSB distance by default is 10
    // 2047 = 0b11111111111, ie. truncation to 11bp
    lastTenStrand1  = lastTenStrand1 & truncator;
    lastTenStrand2  = lastTenStrand2 & truncator;
    lastTenTracked1 = lastTenTracked1 & truncator;
    lastTenTracked2 = lastTenTracked2 & truncator;
    // keep counters to ten binary places
    count1 = count_bytes(lastTenStrand1);
    count2 = count_bytes(lastTenStrand2);

    // Nightmare of if statements
    if(fDamageRecords.at(ii)->fStrand1DirectDmg)
    {
      strand1Direct = true;
      classification->fDirectBreaks++;
    }
    if(fDamageRecords.at(ii)->fStrand1IndirectDmg)
    {
      strand1Indirect = true;
      classification->fIndirectBreaks++;
    }
    if(fDamageRecords.at(ii)->fStrand2DirectDmg)
    {
      strand2Direct = true;
      classification->fDirectBreaks++;
    }
    if(fDamageRecords.at(ii)->fStrand2IndirectDmg)
    {
      strand2Indirect = true;
      classification->fIndirectBreaks++;
    }
    if(fDamageRecords.at(ii)->fbp1DirectDmg)
    {
      // classification->fDirectBreaks++; //???
    }
    if(fDamageRecords.at(ii)->fBp1IndirectDmg)
    {
      if(fDamageRecords.at(ii)->fbp1InducedBreak)
      {
        // Counts as a hit only if it breaks a strand
        classification->fIndirectBreaks++;
        strand1Indirect = true;
        classification->fInducedBreaks++;
      }
    }
    if(fDamageRecords.at(ii)->fbp2DirectDmg)
    {
      // classification->fDirectBreaks++; //???
    }
    if(fDamageRecords.at(ii)->fBp2IndirectDmg)
    {
      if(fDamageRecords.at(ii)->fbp2InducedBreak)
      {
        // Counts as a hit only if it breaks a strand
        classification->fIndirectBreaks++;
        strand2Indirect = true;
        classification->fInducedBreaks++;
      }
    }

    G4bool s1IndirectBreak = fDamageRecords.at(ii)->fStrand1IndirectDmg ||
                             fDamageRecords.at(ii)->fbp1InducedBreak;
    G4bool s2IndirectBreak = fDamageRecords.at(ii)->fStrand2IndirectDmg ||
                             fDamageRecords.at(ii)->fbp2InducedBreak;

    // strand 1 hit
    G4bool strand1hit =
      fDamageRecords.at(ii)->fStrand1DirectDmg || s1IndirectBreak;
    G4bool strand2hit =
      fDamageRecords.at(ii)->fStrand2DirectDmg || s2IndirectBreak;
    G4bool bp1hit = fDamageRecords.at(ii)->fbp1DirectDmg ||
                    fDamageRecords.at(ii)->fBp1IndirectDmg;
    G4bool bp2hit = fDamageRecords.at(ii)->fbp2DirectDmg ||
                    fDamageRecords.at(ii)->fBp2IndirectDmg;
    strandDamage += ((G4int) strand1hit + (G4int) strand2hit);
    baseDamage += ((G4int) bp1hit + (G4int) bp2hit);

    // Update the damage type counters
    // strand 1
    if(s1IndirectBreak)
    {
      lastTenIndirectStrand1++;
    }
    if(fDamageRecords.at(ii)->fStrand1DirectDmg)
    {
      lastTenDirectStrand1++;
    }
    // strand 2
    if(s2IndirectBreak)
    {
      lastTenIndirectStrand2++;
    }
    if(fDamageRecords.at(ii)->fStrand2DirectDmg)
      lastTenDirectStrand2++;

    oldnDSB = nDSB + nDSBPlus;
    if(strand1hit && strand2hit)
    {
      nDSB++;
      if(count1 || count2)
      {
        nDSBPlus++;
      }
      lastTenStrand1++;
      lastTenStrand2++;
    }
    else if(strand1hit)
    {
      if(lastTenTracked2)
      {
        nDSB++;
        lastTenTracked2 =
          (1 << (utility::Fls(lastTenTracked2) - 1)) ^ lastTenTracked2;
        if((count2 > 1) || (count1))
        {
          nDSBPlus++;
        }
      }
      else if(count1 && count2)
      {
        nDSBPlus++;
      }
      else
      {
        nSSB++;
      }
      SSB1 = true;
      lastTenTracked1++;
      lastTenStrand1++;
    }
    else if(strand2hit)
    {
      if(lastTenTracked1)
      {
        nDSB++;
        lastTenTracked1 =
          (1 << (utility::Fls(lastTenTracked1) - 1)) ^ lastTenTracked1;
        if((count1 > 1) || (count2))
        {
          nDSBPlus++;
        }
      }
      else if(count1 && count2)
      {
        nDSBPlus++;
      }
      else
      {
        nSSB++;
      }
      SSB2 = true;
      lastTenStrand2++;
      lastTenTracked2++;
    }
    if(oldnDSB != (nDSB + nDSBPlus))
    {
      // we have had a DSB, time to classify its source.
      // DSB
      if(lastTenDirectStrand1 && lastTenDirectStrand2)
      {
        nDSBd = true;
        if(lastTenIndirectStrand1 || lastTenIndirectStrand2)
          nDSBm = true;
      }
      if(lastTenIndirectStrand1 && lastTenIndirectStrand2)
      {
        nDSBi = true;
        if(lastTenDirectStrand1 || lastTenDirectStrand2)
          nDSBh = true;
        if(lastTenDirectStrand1 && lastTenDirectStrand2)
          nDSBm = true;
      }
      if((lastTenDirectStrand1 && lastTenIndirectStrand2) ||
         (lastTenIndirectStrand1 && lastTenDirectStrand2))
      {
        nDSBh = true;
      }
    }
  }

  // Return based on order of severity

  // G4String complexity = "None";

  complexityEnum complexity = NoneComplexity;
  if(nDSB > 1)
  {
    complexity = DSBplusplus;
  }
  else if(nDSBPlus != 0)
  {
    complexity = DSBplus;
  }
  else if(nDSB != 0)
  {
    complexity = DSB;
  }
  else if((nSSB != 0) && SSB1 && SSB2)
  {
    complexity = twoSSB;
  }
  else if(nSSB > 1)
  {
    complexity = SSBplus;
  }
  else if(nSSB != 0)
  {
    complexity = SSB;
  }
  else
  {
    complexity = NoneComplexity;
  }
  G4bool direct     = (strand1Direct || strand2Direct);
  G4bool indirect   = (strand1Indirect || strand2Indirect);
  sourceEnum source = undefined;
  if(nDSB == 0)
  {
    if(direct && indirect)
    {
      source = SSBm;
    }
    else if(direct)
    {
      source = SSBd;
    }
    else if(indirect)
    {
      source = SSBi;
    }
    else
    {
      source = undefined;
    }
  }
  else
  {
    if(nDSBm != 0)
    {
      source = DSBm;
    }
    else if((nDSBd != 0) && (nDSBi != 0))
    {
      source = DSBm;
    }
    else if((nDSBd != 0) && (nDSBh != 0))
    {
      source = DSBm;
    }
    else if(nDSBd != 0)
    {
      source = DSBd;
    }
    else if(nDSBh != 0)
    {
      source = DSBh;  // Catches DSBi && DSBh
    }
    else if(nDSBi != 0)
    {
      source = DSBi;
    }
    else
    {
      G4ExceptionDescription errmsg;
      this->PrintRecord("", dsbDistance);
      errmsg << "Error in source classification routine" << G4endl;
      errmsg << "I think this is " << complexity << "/"
             << "???" << G4endl;
      G4Exception("DamageRecord::GetClassification", "ERR_BAD_CLASSIFICATION",
                  JustWarning, errmsg);
      source = undefined;
    }
  }

  if(((complexity == NoneComplexity) && (source != undefined)) ||
     ((complexity != NoneComplexity) && (source == undefined)))
  {
    G4ExceptionDescription errmsg;
    errmsg << "You can't have a complexity without a source etc." << G4endl;
    errmsg << "I think this is " << complexity << "/" << source << G4endl;
    this->PrintRecord("", dsbDistance);
    G4Exception("DamageRecord::GetClassification", "ERR_BAD_CLASSIFICATION",
                JustWarning, errmsg);
  }

  classification->fComplexity = complexity;
  classification->fSource     = source;
  classification->fbaseDmg    = baseDamage;
  classification->fStrandDmg  = strandDamage;

  return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// add a contrived test record to the damage record
// params: s1 -> strand 1 damage code
//     b1 -> base 1 damage code
//     b2 -> base 2 damage code
//     s2 -> strand 2 damage code
// codes:
// param < 0: Nothing; param == 0 indirect damage;
// param == 1 direct hit, no damage; param > 1: Direct hit physical damage
void DamageRecord::AddTestDamage(G4int s1, G4int b1, G4int b2, G4int s2)
{
  auto* bp = new BasePairDamageRecord;
  if(s1 == 0)
  {
    bp->fStrand1IndirectDmg = true;
  }
  if(s1 > 0)
  {
    bp->fStrand1Energy = 10 * eV;
  }
  if(s1 > 1)
  {
    bp->fStrand1DirectDmg = true;
  }

  if(b1 == 0)
  {
    bp->fBp1IndirectDmg = true;
  }
  if(b1 > 0)
  {
    bp->fBp1Energy = 10 * eV;
  }
  if(b1 > 1)
  {
    bp->fbp1DirectDmg = true;
  }

  if(s2 == 0)
  {
    bp->fStrand2IndirectDmg = true;
  }
  if(s2 > 0)
  {
    bp->fStrand2Energy = 10 * eV;
  }
  if(s2 > 1)
  {
    bp->fStrand2DirectDmg = true;
  }

  if(b2 == 0)
  {
    bp->fBp2IndirectDmg = true;
  }
  if(b2 > 0)
  {
    bp->fBp2Energy = 10 * eV;
  }
  if(b2 > 1)
  {
    bp->fbp2DirectDmg = true;
  }

  this->AddBasePairDamage(bp, G4ThreeVector(1, 1, G4UniformRand()));
}

// Unit test for classification
void AnalysisManager::TestClassification()
{
  G4cout << G4endl;
  G4cout << "=================================" << G4endl;
  G4cout << "Classification Unit Test Begins" << G4endl;

  // None
  auto* theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(1, -1, -1, 1);
  theRecord->AddEmptyBPDamage(5);
  theRecord->AddTestDamage(-1, -1, -1, 1);
  theRecord->AddEmptyBPDamage(5);
  theRecord->AddTestDamage(-1, -1, -1, 1);
  theRecord->AddEmptyBPDamage(10);
  DamageClassification* classification = theRecord->GetClassification();
  G4cout << "Test None" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == NoneComplexity);
  assert(classification->fSource == sourceEnum::undefined);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // None, base damage
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(1, -1, 0, 1);
  theRecord->AddEmptyBPDamage(5);
  theRecord->AddTestDamage(-1, -1, -1, 1);
  theRecord->AddEmptyBPDamage(5);
  theRecord->AddTestDamage(-1, 0, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test None" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  G4cout << "indirect Breaks = " << classification->fIndirectBreaks << G4endl;
  G4cout << "base Damage = " << classification->fbaseDmg << G4endl;
  assert(classification->fComplexity == NoneComplexity);
  assert(classification->fSource == sourceEnum::undefined);
  assert(classification->fIndirectBreaks == 0);
  assert(classification->fDirectBreaks == 0);
  assert(classification->fbaseDmg == 2);
  assert(classification->fStrandDmg == 0);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBd" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == SSB);
  assert(classification->fSource == SSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBi
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, 0, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBi" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == SSB);
  assert(classification->fSource == SSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBm / SSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, 0, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, 0, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBm / SSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == SSBplus);
  assert(classification->fSource == SSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBm / 2SSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBd / 2SSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == twoSSB);
  assert(classification->fSource == SSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBi / 2SSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBi / 2SSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == twoSSB);
  assert(classification->fSource == SSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBd / SSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, 0, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBd / SSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == SSBplus);
  assert(classification->fSource == SSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBm / SSB2
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBm / 2SSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == twoSSB);
  assert(classification->fSource == SSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBd / SSB2
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBd / 2SSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == twoSSB);
  assert(classification->fSource == SSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // SSBd / SSB2
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(16);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBd / 2SSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == twoSSB);
  assert(classification->fSource == SSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  G4cout << G4endl;

  // Induced Breaks DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);

  auto* bp             = new BasePairDamageRecord();
  bp->fBp1IndirectDmg  = true;
  bp->fbp1InducedBreak = true;
  theRecord->AddBasePairDamage(bp, G4ThreeVector(0, 0, 0));

  bp                   = new BasePairDamageRecord();
  bp->fBp2IndirectDmg  = true;
  bp->fbp2InducedBreak = true;
  theRecord->AddBasePairDamage(bp, G4ThreeVector(0, 0, 0));

  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB induced" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // Induced Breaks SSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);

  bp                   = new BasePairDamageRecord();
  bp->fBp2IndirectDmg  = true;
  bp->fbp2InducedBreak = true;
  theRecord->AddBasePairDamage(bp, G4ThreeVector(0, 0, 0));

  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test SSBi/SSB induced" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == SSB);
  assert(classification->fSource == SSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(9);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi/ DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, 0);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBh / DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBh);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBh / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, 2);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBh);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBm / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, 2);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBh / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBh);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSB);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBh / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBh/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBh);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBm / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBm/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB+
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBd / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(15);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(2, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBd/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBd);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBh / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddTestDamage(0, -1, -1, 0);
  theRecord->AddTestDamage(0, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBh/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBh);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBm / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddTestDamage(2, -1, -1, 0);
  theRecord->AddTestDamage(0, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 2);
  theRecord->AddEmptyBPDamage(12);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddEmptyBPDamage(6);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBm/DSB++" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplusplus);
  assert(classification->fSource == DSBm);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  // DSBi / DSB++
  theRecord = new DamageRecord("Test", 0, 0, 0);
  theRecord->AddEmptyBPDamage(10);
  theRecord->AddTestDamage(-1, 0, -1, -1);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, -1);
  theRecord->AddTestDamage(-1, -1, -1, 0);
  theRecord->AddTestDamage(-1, -1, -1, -1);
  theRecord->AddTestDamage(-1, -1, -1, -1);
  theRecord->AddTestDamage(0, -1, -1, -1);
  theRecord->AddEmptyBPDamage(10);
  classification = theRecord->GetClassification();
  G4cout << "Test DSBi/DSB+" << G4endl;
  theRecord->PrintRecord("", 10);
  G4cout << "Classification: " << classification->fSource << "/"
         << classification->fComplexity << G4endl;
  assert(classification->fComplexity == DSBplus);
  assert(classification->fSource == DSBi);
  delete classification;
  delete theRecord;
  G4cout << "---------------------------------" << G4endl;

  G4cout << "Classification Unit Test Complete" << G4endl;
  G4cout << "=================================" << G4endl;

  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageRecord::AddEmptyBPDamage(const int64_t& ii)  // dousatsu
{
  auto basePairNumber = ii;
  while(basePairNumber > 0)
  {
    fDamageRecords.push_back(new BasePairDamageRecord);
    basePairNumber--;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageRecord::AddStrandHit(const G4MoleculeDefinition* mol)
{
  if(mol == G4OH::Definition())
  {
    fOHStrand++;
  }
  else if(mol == G4Electron_aq::Definition())
  {
    fEaqStrand++;
  }
  else if(mol == G4Hydrogen::Definition())
  {
    fHStrand++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageRecord::AddBaseHit(const G4MoleculeDefinition* mol)
{
  // G4cout << "Hit by " << mol->GetParticleName() << G4endl;
  if(mol == fOH)
  {
    fOHBase++;
  }
  else if(mol == fe_aq)
  {
    fEaqBase++;
  }
  else if(mol == fH)
  {
    fHBase++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector DamageRecord::GetMeanPosition() const
{
  G4double x(0.);
  G4double y(0.);
  G4double z(0.);
  for(const auto& fPosition : fPositions)
  {
    x += fPosition.getX();
    y += fPosition.getY();
    z += fPosition.getZ();
  }

  if(fPositions.empty())
  {
    G4ExceptionDescription errmsg;
    errmsg << "fPositions is emply ";
    G4Exception("DamageRecord", "ERR_INVALID_PROB", FatalException, errmsg);
    return G4ThreeVector(0., 0., 0.);
  }
  else
  {
    G4double factor = 1.0 / fPositions.size();
    return G4ThreeVector(x * factor, y * factor, z * factor);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DamageRecord::GetMeanDistance() const
{
  G4double d = 0.;
  for(auto it = fPositions.begin(); it != fPositions.end(); ++it)
  {
    for(auto jt = it + 1; jt != fPositions.end(); ++jt)
    {
      d += ((*jt) - (*it)).mag();
    }
  }
  G4int number = fPositions.size();
  return d / number;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DamageRecord::GetEnergy() const
{
  G4double en = 0;
  for(auto bp : fDamageRecords)
  {
    en = en + bp->fBp1Energy + bp->fBp2Energy + bp->fStrand1Energy +
         bp->fStrand2Energy;
  }
  return en;
}