{
//
// ROOT MACRO FOR TESTING G4/HETC++ DATA
//
// read data from an ascii file and
// create a root file with an histogram and an ntuple.

gROOT->Reset();
#include <iostream.h>

ifstream in;

in.open("gaurn.out", ios::in); // data must be in current directory

Float_t x;
Int_t nlines = 0;

gBenchmark->Start("fit1");
TH1F *h1 = new TH1F("h1","gaurn",150,-4,4);
TF1 *func = new TF1("func","gaus(0)",-4,4);
func->SetParameters(10,4,1,20);

TCanvas* c1=new TCanvas("c1", "new canvas");
c1->SetFillColor(0);
while (1) {
  in >> x;
  if (!in.good()) break;
  h1->Fill(x);
  nlines++;
}
in.close();
printf(" found %d points\n",nlines);
h1->Fit("func");
h1->GetXaxis()->SetTitle("return value");
h1->GetXaxis()->CenterTitle();
h1->GetYaxis()->SetTitle("number of hits");
h1->GetYaxis()->CenterTitle();

gBenchmark->Show("fit1");


//ntuple->Draw("x"); // show results


} // end of macro

