{
   
  gROOT->Reset();
  gROOT->SetStyle("clearRetro");

#include "Riostream.h"
  float MeV=0.001;
  ifstream in;
  in.open("../data/b2.out");

  int runId, eventId, particleId, fragmentA, fragmentZ;
  float kineticEnergy, momX, momY, momZ, exitationEnergy;
  int nlines = 0;
  TFile *f = new TFile("analyzeEvents.root","RECREATE");
  TH1F *hpE = new TH1F("hpE","proton kinetic energy",1000,0,100);
  TH1F *hnE = new TH1F("hnE","neutron kinetic energy",1000,0,100);
  TNtuple *ntuple = new TNtuple("ntuple","foorified data from cascade.cc  ouput","runId:eventId:particleId:kineticEnergy:momX:momY:momZ:fragmentA:fragmentZ:exitationEnergy");

  while (1) {
    in >> runId >> eventId >> particleId >> kineticEnergy >> momX >> momY >> momZ >> fragmentA >> fragmentZ >> exitationEnergy;
    if (!in.good()) break;
    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",runId, eventId, particleId);
    if (particleId==1) hpE->Fill(kineticEnergy/MeV);
    if (particleId==0) hnE->Fill(kineticEnergy/MeV);
    ntuple->Fill(runId, eventId, particleId, kineticEnergy, momX, momY, momZ, fragmentA, fragmentZ, exitationEnergy);
    nlines++;
  }
  printf(" found %d lines \n",nlines);

  in.close();

  f->Write();
}
