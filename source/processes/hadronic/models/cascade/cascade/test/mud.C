{
//
// ROOT MACRO FOR TESTING G4/HETC++ DATA
//
// read data from an ascii file and
// create a root file with an histogram and an ntuple.

gROOT->Reset();

#include <iostream.h>

ifstream in;

in.open("mud.out", ios::in); // data must be in current directory

Float_t x,y,z;
Int_t nlines = 0;

TFile *f = new TFile("mud.root","RECREATE");
TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y:z");

while (1) {
  in >> x >> y >> z;
  if (!in.good()) break;
  if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
  h1->Fill(x);
  ntuple->Fill(x,y,z);
  nlines++;
}
printf(" found %d points\n",nlines);

in.close();

f->Write();


ntuple->Draw("x"); // show results
} // end of macro


