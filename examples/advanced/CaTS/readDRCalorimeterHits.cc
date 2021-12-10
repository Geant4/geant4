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
/// \file readDRCalorimeterHits.cc
/// \brief example how to read the  CaTS::DRCalorimeterHits
//
// Root headers
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
// project headers
#include "Event.hh"
#include "DRCalorimeterHit.hh"

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
  TFile fo(argv[1]);
  fo.GetListOfKeys()->Print();
  Event* event = new Event();
  TTree* Tevt  = (TTree*) fo.Get("Events");
  Tevt->SetBranchAddress("event.", &event);
  TBranch* fevtbranch = Tevt->GetBranch("event.");
  Int_t nevent        = fevtbranch->GetEntries();
  G4cout << " Nr. of Events:  " << nevent << G4endl;
  double max                 = 0.;
  double min                 = 1000000;
  double nmax                = 0.;
  double nmin                = 1000000000;
  std::string CollectionName = argv[3];
  CollectionName             = CollectionName + "_DRCalorimeter_HC";
  for(Int_t i = 0; i < nevent; i++)
  {
    fevtbranch->GetEntry(i);
    auto* hcmap = event->GetHCMap();
    for(const auto& ele : *hcmap)
    {
      if(ele.first.compare(CollectionName) == 0)
      {
        auto hits    = ele.second;
        G4int NbHits = hits.size();
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          DRCalorimeterHit* drcaloHit =
            dynamic_cast<DRCalorimeterHit*>(hits.at(ii));
          const double ed = drcaloHit->GetEdep();
          if(ed > max)
            max = ed;
          if(ed < min)
            min = ed;
          const unsigned int nceren = drcaloHit->GetNceren();
          if(nceren > nmax)
            nmax = nceren;
          if(nceren < nmin)
            nmin = nceren;
        }
      }
    }
  }
  outfile->cd();
  TH1F* hedep   = new TH1F("energy", "edep", 100, min, max);
  TH1F* hnceren = new TH1F("nceren", "nceren", 100, nmin, nmax);
  for(Int_t i = 0; i < nevent; i++)
  {
    fevtbranch->GetEntry(i);
    auto* hcmap = event->GetHCMap();
    for(const auto& ele : *hcmap)
    {
      if(ele.first.compare(CollectionName) == 0)
      {
        auto hits    = ele.second;
        G4int NbHits = hits.size();
        for(G4int ii = 0; ii < NbHits; ii++)
        {
          DRCalorimeterHit* drcaloHit =
            dynamic_cast<DRCalorimeterHit*>(hits.at(ii));
          hedep->Fill(drcaloHit->GetEdep());
          hnceren->Fill(drcaloHit->GetNceren());
        }
      }
    }
  }
  // hedep->Fit("gaus");
  // hnceren->Fit("gaus");
  outfile->Write();
}
