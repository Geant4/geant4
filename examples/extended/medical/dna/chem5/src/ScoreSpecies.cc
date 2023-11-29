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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// Phys. Med. Biol. 63(10) (2018) 105014-12pp
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// ScoreSpecies.cc
//
#include "ScoreSpecies.hh"

#include "G4UnitsTable.hh"
#include <G4MolecularConfiguration.hh>
#include <G4MoleculeCounter.hh>
#include "G4Event.hh"
#include <G4SystemOfUnits.hh>
#include <globals.hh>
#include <G4EventManager.hh>
#include <iomanip>

/**
 \file ScoreSpecies.cc
 \class ScoreSpecies
  This is a primitive scorer class for molecular species.
  The number of species is recorded for all times (log spaced). It 
  also scores the energy deposition in order to compute the 
  radiochemical yields.
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScoreSpecies::ScoreSpecies(G4String name, G4int depth)
: G4VPrimitiveScorer(name,depth),
  fEdep(0),
  fHCID(-1),
  fEvtMap(0)
{
  fNEvent = 0;
  G4double tMin = 1.0 * CLHEP::picosecond;
  G4double tMax = 999999 * CLHEP::picosecond;
  G4double tLogMin = std::log10(tMin);
  G4double tLogMax = std::log10(tMax);
  G4int tBins = 50;
  for ( int i = 0; i <= tBins; i++ )  
    AddTimeToRecord(std::pow(10., tLogMin + i*(tLogMax-tLogMin)/tBins));
  
  fEdep = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScoreSpecies::~ScoreSpecies()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ScoreSpecies::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ( edep == 0. ) return FALSE;

  edep *= aStep->GetPreStepPoint()->GetWeight();
  G4int  index = GetIndex(aStep);
  fEvtMap->add(index,edep);
  fEdep+=edep;

  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::Initialize(G4HCofThisEvent* HCE)
{
  fEvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                    GetName());

  if(fHCID < 0)
  {
    fHCID = GetCollectionID(0);
  }

  HCE->AddHitsCollection(fHCID, (G4VHitsCollection*)fEvtMap);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::EndOfEvent(G4HCofThisEvent*)
{
  if(G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted())
  {
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }
  
  auto species = G4MoleculeCounter::Instance()->GetRecordedMolecules();

  if(species.get() == 0 || species->size() == 0)
  {
    G4cout << "No molecule recorded, energy deposited= "
           << G4BestUnit(fEdep, "Energy") << G4endl;
    ++fNEvent;
    fEdep = 0.;
    G4MoleculeCounter::Instance()->ResetCounter();
    return;
  }
  
  for(auto molecule: *species)
  {
    for(auto time_mol: fTimeToRecord)
    {
      double n_mol =
          G4MoleculeCounter::Instance()->GetNMoleculesAtTime(molecule,
                                                             time_mol);

      if(n_mol < 0)
      {
        G4cerr << "N molecules not valid < 0 " << G4endl;
        G4Exception("","N<0",FatalException,"");
      }

      SpeciesInfo& molInfo = fSpeciesInfoPerTime[time_mol][molecule];
      molInfo.fNumber += n_mol;
      double gValue = (n_mol/(fEdep/eV)) * 100.;
      molInfo.fG += gValue;
      molInfo.fG2 += gValue*gValue;
    }
  }
  ++fNEvent;
  
  fEdep = 0.;
  G4MoleculeCounter::Instance()->ResetCounter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
ScoreSpecies::AbsorbResultsFromWorkerScorer(G4VPrimitiveScorer* workerScorer)
{
  ScoreSpecies* right =
  dynamic_cast<ScoreSpecies*>(dynamic_cast<G4VPrimitiveScorer*>(workerScorer));

  if(right == 0)
  {
    return;
  }
  if(right == this)
  {
    return;
  }

  SpeciesMap::iterator it_map1 = right->fSpeciesInfoPerTime.begin();
  SpeciesMap::iterator end_map1 = right->fSpeciesInfoPerTime.end();

  for(; it_map1 != end_map1; ++it_map1)
  {
    InnerSpeciesMap& map2 = it_map1->second;
    InnerSpeciesMap::iterator it_map2 = map2.begin();
    InnerSpeciesMap::iterator end_map2 = map2.end();

    for(; it_map2 != end_map2; ++it_map2)
    {
      SpeciesInfo& molInfo =
      fSpeciesInfoPerTime[it_map1->first][it_map2->first] ;
      molInfo.fNumber  += it_map2->second.fNumber;
      molInfo.fG += it_map2->second.fG;
      molInfo.fG2 += it_map2->second.fG2;
    }
  }
  
  fNEvent += right->fNEvent;
  right->fNEvent = 0;
  right->fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::DrawAll()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::PrintAll()
{
   G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
   G4cout << " PrimitiveScorer " << GetName() << G4endl;
   G4cout << " Number of events " << fNEvent << G4endl;
   G4cout << " Number of energy deposition recorded "
          << fEvtMap->entries() << G4endl;
  
  for(auto itr : *fEvtMap->GetMap()) {
     G4cout << "  copy no.: " << itr.first
     << "  energy deposit: "
     << *(itr.second)/GetUnitValue()
     << " [" << GetUnit()<<"]"
     << G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::ASCII()
{
  std::ofstream out("Species.txt");
  if(!out) return;

  out << "# Time [ps]  G-value (/100 eV)  RMS   Molecule" << G4endl;

  std::map<G4String, std::map<G4double, std::pair<G4double,G4double>>> mol;
 
  for(auto it_map1: fSpeciesInfoPerTime)
  {
    InnerSpeciesMap& map2 = it_map1.second;
    G4double time = it_map1.first/ps;
    for(auto it_map2: map2)
    {
      G4double G = it_map2.second.fG;
      G4double G2 = it_map2.second.fG2;
      G4double N = fNEvent;
      G /= N;
      G2 = std::sqrt( N/(N-1) * ( G2/N - G*G) );
      mol[it_map2.first->GetName()][time]=std::make_pair(G,G2);
    }
  }

  for ( auto it1 : mol )
    for ( auto it2 : it1.second )
      out << std::setw(12) << it2.first << std::setw(12) << it2.second.first 
          << std::setw(12) << it2.second.second << std::setw(12) 
          << std::setw(12) << it1.first << G4endl; 

  out.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScoreSpecies::OutputAndClear()
{
  if(G4Threading::IsWorkerThread()) return;

  //----------------------------------------------------------------------------
  // Save results in ASCII format
  ASCII();

  fNEvent = 0;
  fSpeciesInfoPerTime.clear();
}

