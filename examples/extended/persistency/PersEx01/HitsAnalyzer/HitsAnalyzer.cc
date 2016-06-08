//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include <stdlib.h>
#include <getopt.h>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TRandom.h"

#include "Pers01CalorHit.hh"
//#include "Pers01CalorHitRoot.hh"
#include "Pers01CalorHitRootIO.hh"
#include "G4PersistencyCenter.hh"

//////////////////////////////////////////////////////////////////////

void usage(char*);

//////////////////////////////////////////////////////////////////////

void usage(char* myname)
{
//usage: HitsAnalyzer [-h] [-v] [-r|w] 
//                    [-f <file>] [-o <histofile>] [-n <n>] [-s <n>]

  std::cerr << std::endl
            << "usage: " << myname
            << " [-h] [-v] [-r|w]" << std::endl
            << "                     "
            << " [-f <file>] [-o <histofile>]"
            << " [-n <n>] [-s <n>]" << std::endl
            << "   -h: print this menu"  << std::endl
            << "   -v: verbose output"   << std::endl
            << "   -r: read the file"    << std::endl
            << "   -w: write the file"   << std::endl
            << "   -f: file name"        << std::endl
            << "   -n: number of events" << std::endl
            << "   -s: split level"      << std::endl
            << "   -N: Total nodes (Gfarm option)" << std::endl
            << "   -I: Node index (Gfarm option)" << std::endl
            << std::endl;
  exit(1);
}

int main(int argc, char **argv)
{
  Int_t nevent = 15;  
  Int_t split  = 1;       // by default, split Event in sub branches
  Int_t write  = 1;       // by default the event is filled
  Int_t read   = 0;
  Int_t verbose = 0;

  std::string filename = "Event.root";
  std::string histofile = "Histo.root";

  int c;
  while ((c = ::getopt(argc, argv, "hvrwf:n:s:")) != EOF) {
    switch (c) {
      case 'v':
        verbose = 1;
        break;
      case 'r':
        read  = 1;
        write = 0;
        break;
      case 'w':
        read  = 0;
        write = 1;
        break;
      case 'f':
        filename = optarg;
        break;
      case 'o':
        histofile = optarg;
        break;
      case 'n':
        nevent = (int) ::atoi(optarg);
        break;
      case 's':
        split = (int) ::atoi(optarg);
        break;
      case 'h':
      default:
        usage(argv[0]);
        break;
    }
  }

  if (verbose) {
    if (read)  cout << "*** Input File: " << filename << endl;
    if (write) cout << "*** Output File: " << filename << endl;
    cout << "*** Histogram File: " << histofile << endl;
  }

  Int_t ev;
  Int_t printev = 10000;

  cout << endl << " *** Testing Pers01CalorHitRootIO *** " << endl << endl;

  string detName = "CalorSD";
  string colName = "CalCollection";
  string obj     = "Hits";
  string pmName  = "ROOT";

  G4PersistencyCenter* pc = G4PersistencyCenter::GetPersistencyCenter();
  G4PersistencyManager* pm = pc->GetPersistencyManager(pmName);
  pm = pm->Create();
  pc->SetPersistencyManager(pm,pmName);
  pc->SetReadFile("Hits",filename);
  pc->SetWriteFile("Hits",filename);

  Pers01CalorHitRootIO* ioman = new Pers01CalorHitRootIO(detName,colName);
  if (verbose) ioman->SetVerboseLevel(5);

  Pers01CalorHitsCollection* hc = 0;

//         Read case
  if (read) {
    std::cout << "Reading ..." << std::endl;

    pm->TransactionManager()->StartRead();

    TFile* hf = new TFile(histofile.c_str(), "RECREATE");
    TObjArray Hlist(0);
    TH1F* h1 = new TH1F("hist1", "Pers01CalorHit EdepAbs", 100, 0, 4.4);
    TH1F* h2 = new TH1F("hist2", "Pers01CalorHit TrackLengthAbs", 100, 0, 4.4);
    TH1F* h3 = new TH1F("hist3", "Pers01CalorHit EdepGap", 100, 0, 4.4);
    TH1F* h4 = new TH1F("hist4", "Pers01CalorHit TrackLengthGap", 100, 0, 4.4);
    Hlist.Add(h1);
    Hlist.Add(h2);
    Hlist.Add(h3);
    Hlist.Add(h4);

    G4VHitsCollection* ahc;

    string file = pc->CurrentReadFile(obj);
    if ( pm->TransactionManager()->SelectReadFile(obj, file) ) {

      ev = 0;
      while ( ioman->Retrieve(ahc) ) {
        hc = (Pers01CalorHitsCollection*) ahc;
        assert(ahc);

        if (verbose || ev%printev == 0) printf("event:%d\n",ev);
        if ( ev++ > nevent ) break;
    
        // analyze transient Pers01CalorHit
        hc = (Pers01CalorHitsCollection*) ahc;
        int n = hc->entries();
        for ( int i = 0; i < n; i++ ) {
    
          Pers01CalorHit* h = (*hc)[i];
          if ( h != 0 ) {
            h1->Fill(h->GetEdepAbs());
            h2->Fill(h->GetTrakAbs());
            h3->Fill(h->GetEdepGap());
            h4->Fill(h->GetTrakGap());
          }
        }
        delete ahc;
      }
    }

    pm->TransactionManager()->Commit();
    hf->Write();

  } else {
//         Write case
    std::cout << "Writing ..." << std::endl;
  
    pm->TransactionManager()->StartUpdate();

    for (ev = 0; ev < nevent; ev++) {
      if (verbose || ev%printev == 0) printf("event:%d\n",ev);

      Float_t edep = 600*gRandom->Rndm(1);
  
      Pers01CalorHit* hit = new Pers01CalorHit();
      hit->AddAbs(edep,2*edep);
      hit->AddGap(edep,2*edep);

      if (verbose) cout << " EdepAbs = " << hit->GetEdepAbs() << endl;

      hc = new Pers01CalorHitsCollection(detName,colName);
      hc->insert(hit);

      string file = pc->CurrentWriteFile(obj);
      if ( pm->TransactionManager()->SelectWriteFile(obj, file) ) {
        ioman->Store(hc);
      }
      delete hc;
    }

    pm->TransactionManager()->Commit();
  }

  cout << "Done." << endl;

  // delete hits collection, IO manager
  // delete pc;
  delete ioman;

}

