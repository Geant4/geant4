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
// CaTS (Calorimetry and Tracking Simulation)
//
// Authors: Hans Wenzel and Soon Yung Jun
//          (Fermi National Accelerator Laboratory)
//
// History: October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file RootIO.hh
/// \brief Definition of the CaTS::RootIO class

#ifdef WITH_ROOT
#pragma once

#include <G4Types.hh>
#include <G4String.hh>
#include "G4ThreadLocalSingleton.hh"

#include "G4VHit.hh"
#include "G4VHitsCollection.hh"
#include "Event.hh"

class G4Event;
class TBranch;
class TFile;
class TTree;

class RootIO
{
  friend class G4ThreadLocalSingleton<RootIO>;

 public:
  static RootIO* GetInstance();

  ~RootIO() = default;
  void Write(const G4Event* event);
  void Merge();
  void Close();

 private:
  // Private constructor used only for the static instance
  RootIO();

  // Split a string using the delimiter character
  std::vector<G4String> Split(const G4String& s, char delim);

  // Add hits from sensitive detectors into the CaTS event data
  template<typename T>
  void AddHits(G4VHitsCollection* hc, Event* event_data);

private:
  G4String fFileName;
  TBranch* fEvtBranch{ nullptr };
  TFile* fFile{ nullptr };
  TTree* fTree{ nullptr };
};

template<typename T>
void RootIO::AddHits(G4VHitsCollection* hc, Event* event_data)
{
  std::map<G4String, std::vector< G4VHit*>> *hcmap = event_data->GetHCMap();
  std::vector<G4VHit*> hitsVector;
  G4int NbHits = hc->GetSize();
  G4String hcname = hc->GetName();

  for (G4int ii = 0; ii < NbHits; ++ii)
  {
    G4VHit* hit = hc->GetHit(ii);
    T* tpcHit = dynamic_cast<T*> (hit);
    hitsVector.push_back(tpcHit);
  }
  hcmap->insert(std::make_pair(hcname, hitsVector));
}

#endif /* WITH_ROOT */
