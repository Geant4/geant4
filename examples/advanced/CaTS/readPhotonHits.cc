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
//---------------------------------------------------------------------
//*            |\___/|                                                *
//*            )     (                                                *
//*           =\     /=                                               *
//*             )===(                                                 *
//*            /     \     CaTS: Calorimeter and Tracker Simulation   *
//*            |     |     is a flexible and extend-able framework    *
//*           /       \    for the simulation of various detector     *
//*	      \       /    systems                                    *
//*            \__  _/     https://github.com/hanswenzel/CaTS         *
//*	         ( (                                                  *
//*	          ) )                                                 *
//*              (_(                                                  *
//* CaTS also serves as an example that demonstrates how to use       *
//* opticks from within Geant4 for the creation and propagation of    *
//* optical photons.                                                  *
//* see https://bitbucket.org/simoncblyth/opticks.git).               *
//* Ascii Art by Joan Stark: https://www.asciiworld.com/-Cats-2-.html *
//---------------------------------------------------------------------
//
/// \file readPhotonHits.cc
/// \brief example how to read the  CaTS::PhotonHits
//
// Root headers
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
// Project headers
#include "Event.hh"
#include "lArTPCHit.hh"
#include "PhotonHit.hh"

int main(int argc, char** argv)
{
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCaTSClassesDict");
  if(argc < 4)
  {
    G4cout << "Program requires 3 arguments: name of input file, name of "
              "output file, Volume that sensitive detector is attached to"
           << G4endl;
    exit(1);
  }
  TFile* outfile = new TFile(argv[2], "RECREATE");
  outfile->cd();
  TH2F* pos2  = new TH2F("position", "position of Photon Hits", 400, -1000.,
                        1000., 400, -500, 500);
  TH1F* time  = new TH1F("time", "timing of photon hits", 1000, 0., 250.);
  TH1F* time0 = new TH1F("time0", "timing of photon hits", 1000, 0., 250.);
  TH1F* time1 = new TH1F("time1", "timing of photon hits", 1000, 0., 250.);
  TH1F* time2 = new TH1F("time2", "timing of photon hits", 1000, 0., 250.);
  TH1F* wl = new TH1F("wl", "wavelength of detected photons", 1000, 0., 250.);
  TH1F* np = new TH1F("np", "number of detected photons", 100, 0., 50.);
  TFile fo(argv[1]);
  fo.GetListOfKeys()->Print();
  Event* event = new Event();
  TTree* Tevt  = (TTree*) fo.Get("Events");
  Tevt->SetBranchAddress("event.", &event);
  TBranch* fevtbranch = Tevt->GetBranch("event.");
  Int_t nevent        = fevtbranch->GetEntries();
  G4cout << "Nr. of Events:  " << nevent << G4endl;
  std::string CollectionName = argv[3];
  CollectionName             = CollectionName + "_Photondetector_HC";
  for(Int_t i = 0; i < nevent; i++)
  {
    fevtbranch->GetEntry(i);
    auto* hcmap = event->GetHCMap();
    for(const auto& ele : *hcmap)
    {
      auto hits = ele.second;
      if(ele.first.compare(CollectionName) == 0)
      {
        auto hits    = ele.second;
        G4int NbHits = hits.size();
        G4cout << "Event: " << i << "  Number of Hits:  " << NbHits << G4endl;
        np->Fill(NbHits);
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          PhotonHit* photonHit = dynamic_cast<PhotonHit*>(hits.at(ii));
          time->Fill(photonHit->GetTime());
          wl->Fill(photonHit->GetWavelength());
          if(photonHit->GetPosition().getZ() < -100.)
            time0->Fill(photonHit->GetTime());
          if(photonHit->GetPosition().getZ() > -100. &&
             photonHit->GetPosition().getZ() < 100)
            time1->Fill(photonHit->GetTime());
          if(photonHit->GetPosition().getZ() < 100.)
            time2->Fill(photonHit->GetTime());
          pos2->Fill(photonHit->GetPosition().getZ(),
                     photonHit->GetPosition().getY());
        }
      }
    }
  }
  outfile->cd();
  outfile->Write();
}
