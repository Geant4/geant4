{
// Batch analysis for g4 cascade simulation: root -b -q doubleDifferential.C
// Make double differential plots.

// Gray shades 10...19, 0 (from flack to white) 40-49 shades of red
// enum particleType { fragment = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

gROOT->Reset();

printf(" ::: Reading cascade data ...\n");

ifstream in;
in.open("cascade.out", ios::in);

Float_t nEve, typePart, eKin, momX, momY, momZ, nucEx;
Int_t nucA, nucZ,sumB,sumE;
Int_t nlines = 0;

TFile *f = new TFile("cascade.root","RECREATE");
TH1F *h1 = new TH1F("h1","particle types",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from cascade.out","nEve:typePart:eKin:momX:momY:momZ:nucA:nucZ:nucEx");

while (1) {
  in >> nEve >> typePart >> eKin >> momX >> momY >> momZ >> nucA >> nucZ >> nucEx;
  if (!in.good()) break;

  // Show first lines of data
  if (nlines < 2) printf("nEve = %1i, typePart = %1i, eKin = %4f\n",nEve, typePart, eKin);

  // h1->Fill(typePart);
  ntuple->Fill(nEve, typePart, eKin, momX, momY, momZ, nucA, nucZ, nucEx);
  nlines++;
}

printf(" ::: Found %d points\n", nlines);

in.close();
f->Write();

// r(results)
r = new TCanvas("c1","Analysis from cascade.cc ", 200, 10, 700, 780);
r->SetFillColor(0);
r->SetBorderMode(0);

Int_t n = 2; // rows
Int_t m = 3; // columns 
r->Divide(n, m);

Int_t i = 1;
Int_t nm = m * n +1;

while (i < nm) {
  r->cd(i);
  gPad->SetBorderMode(0);
  gPad->Draw();
   gPad->SetLogy();
  //gPad->SetLogx();
  i++;
};

// Change default style 
gStyle->SetStatW(0.30);
gStyle->SetStatH(0.20);
gStyle->SetStatColor(0);
gStyle->SetOptStat(10);

// Proton energy
r->cd(1);
ntuple->SetLineColor(1);
ntuple->SetFillColor(0);
ntuple->Draw("eKin", "typePart==1 ");

r->cd(2);
//ntuple->SetLineColor(1);
//ntuple->SetFillColor(0);
//ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY) > 0.0");
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY) > 0.1");
r->cd(3);
ntuple->SetLineWidth(1);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 0 && sqrt(momX*momX+momY*momY)/momZ < 1");
ntuple->SetLineWidth(2);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 1 && sqrt(momX*momX+momY*momY)/momZ < 2", "same");
ntuple->SetLineWidth(3);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 2 && sqrt(momX*momX+momY*momY)/momZ < 3", "same");
ntuple->SetLineWidth(1);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 3 && sqrt(momX*momX+momY*momY)/momZ < 4", "same");

r->cd(4);
ntuple->SetLineWidth(1); // .18 = 10 deg
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 0.18 && sqrt(momX*momX+momY*momY)/momZ < 0.36");  // 10-20 deg tan(20)=0.36=pt/pl
ntuple->SetLineWidth(2);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 1.19 && sqrt(momX*momX+momY*momY)/momZ < 2.75", "same");  // 50-70

ntuple->SetLineWidth(3);
ntuple->Draw("eKin", "typePart==1 && sqrt(momX*momX+momY*momY)/momZ > 5.67", "same");  // 80-

r->Update();

r->Print("doubleDifferential.eps");

printf(" ::: Done.\n");
}

