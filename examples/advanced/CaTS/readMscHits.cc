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
/// \file readMscHits.cc
/// \brief example how to read the  CaTS::MscHits
//
// Root headers
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1.h"
#include <string>
// project headers
#include "Event.hh"
#include "MscHit.hh"
#include "G4ThreeVector.hh"

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
  TFile fo(argv[1]);
  Event* event = new Event();
  TTree* Tevt  = (TTree*) fo.Get("Events");
  Tevt->SetBranchAddress("event.", &event);
  TBranch* fevtbranch = Tevt->GetBranch("event.");
  Int_t nevent        = fevtbranch->GetEntries();
  outfile->cd();
  TH1F* theta                = new TH1F("angle", "angle", 100, 0.0, 0.15);
  std::string CollectionName = argv[3];
  CollectionName             = CollectionName + "_Msc_HC";
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
          MscHit* mscHit = dynamic_cast<MscHit*>(hits.at(ii));
          G4ThreeVector invec(0, 0, 1.);
          std::cout << std::setprecision(12) << mscHit->GetKinE() << " "
                    << std::setprecision(12) << mscHit->GetMomentum().getX()
                    << " " << std::setprecision(12)
                    << mscHit->GetMomentum().getY() << " "
                    << std::setprecision(12) << mscHit->GetMomentum().getZ()
                    << std::endl;
          theta->Fill(invec.angle(mscHit->GetMomentum()) * 57.3);
        }
      }
    }
  }
  outfile->Write();
}
