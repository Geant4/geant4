{
//simple macro for drawing 1-D histograms
gROOT->Reset();

#include <iostream.h>

ifstream in;

in.open("hist.out", ios::in); // data must be in current directory

Float_t x;
TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x");

while (1) {
  in >> x;
  if (!in.good()) break;
  ntuple->Fill(x,x,x);
}
in.close();


TCanvas* c1=new TCanvas("c1", "new canvas");
c1->SetFillColor(0);
ntuple->Draw("x"); // show results
} // end of macro

