{
// automated analysis for g4 cascade simulation
// root -b -q cascade.C

#include <iostream.h>

printf("::: Reading cascade data ...\n");

gROOT->Reset();

ifstream in;
// we assume a file basic.dat in the current directory
// this file has 3 columns of float data
in.open("cascade.out", ios::in);


Float_t   nEvents, nPart, nProton, nNeutron, nucKinE, pKinE, nKinE, nPi, piKinE;
Int_t nlines = 0;
TFile *f = new TFile("cascade.root","RECREATE");
TH1F *h1 = new TH1F("h1","nucKinE",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from cascade.out","nEvents:nPart:nProton:nNeutron:nucKinE:pKinE:nKinE:nPi:piKinE");

while (1) {
  in >> nEvents >> nPart >> nProton >> nNeutron >> nucKinE >> pKinE >> nKinE >> nPi >> piKinE;
  if (!in.good()) break;
  if (nlines < 5) printf("nPart=%8f, nProton=%8f, nNeutron=%8f\n",nPart,nProton,nNeutron);
  h1->Fill(nucKinE);
  ntuple->Fill(nEvents, nPart, nProton, nNeutron, nucKinE, pKinE, nKinE, nPi, piKinE);
  nlines++;
}
printf(" found %d points\n",nlines);

in.close();

f->Write();

printf("::: Writing ps files ...\n");

TPostScript* ps = new TPostScript("cascade.ps",-111);
   ps->Off();
   ps->On();
   ps->NewPage();
   c1->Update();
   ps->Off();
   ps->Close();
}
