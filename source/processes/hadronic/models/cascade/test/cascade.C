{
// Batch analysis for g4 cascade simulation: root -b -q cascade.C
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
  //  gPad->SetLogy();
  //gPad->SetLogx();
  i++;
}

// Change default style 
gStyle->SetStatW(0.30);
gStyle->SetStatH(0.20);
gStyle->SetStatColor(0);
gStyle->SetOptStat(10);

// Proton energy
r->cd(1);
ntuple->SetLineColor(1);
ntuple->SetFillColor(0);
ntuple->Draw("eKin", "typePart==1");

SHOWEXP = 1;  // Flag to put experimental data to plot.
if (SHOWEXP) {
  // experimental data (run with 3210 events) sigma_rec = 321 mb, original data c= 1
  const Int_t n = 6;
  const Int_t c = 10;
  Double_t exp[n][2] = 
  {{0.003, 20.0 * c},
   {0.005, 14.0 * c},
   { 0.01,  8.0* c}, 
   { 0.02,  5.5* c},
   { 0.04,  4.0* c},
   {0.055,  2.0* c}};

  Double_t e[n], exp[n], measurement[n];

  for (Int_t i = 0; i < n; i++) {
    e[i] = exp[i][0];
    measurement[i] = exp[i][1];
    printf("  %i energy = %f cross-section %f \n", i, e[i], measurement[i]);
  }

  gr0 = new TGraph(n, e, measurement);
  gr0->SetLineColor(1);
  gr0->SetLineWidth(1);
  gr0->SetMarkerColor(1);
  gr0->SetMarkerStyle(21);
  gr0->SetMarkerSize(0.5);
  gr0->SetTitle("");
  gr0->Draw("p");


  TLatex l;
  l->DrawLatex(0.01,  200, "C^{12}_{6}(p,xp) at 62 MeV");
  l->DrawLatex(5,  18, "");
};

// Neutron energy
r->cd(2);

gPad->GetFrame()->SetFillColor(3);

ntuple->SetFillColor(18);
ntuple->Draw("eKin", "typePart==2");

// Pions
r->cd(3);

gPad->GetFrame()->SetFillColor(0);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 

// pi-
ntuple->SetLineStyle(1);
ntuple->SetLineColor(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("eKin", "typePart==3");

// pi+ 
ntuple->SetLineColor(16);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(2); 
//ntuple->Draw("eKin", "typePart==3","same");

// pi0
//ntuple->SetLineStyle(4); // dots
ntuple->SetLineColor(13);
ntuple->SetLineWidth(3); 
//ntuple->Draw("eKin", "typePart==7", "same");

// Fotons
r->cd(4);
gPad->GetFrame()->SetFillColor(0);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 
ntuple->SetLineColor(1);
ntuple->SetFillColor(0);
ntuple->Draw("eKin", "typePart==7");

// Nuclei
r->cd(5);
//gPad->GetFrame()->SetFillColor(0);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 

ntuple->SetLineColor(16);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("nucA", "typePart==0");
ntuple->SetLineColor(13);
ntuple->SetLineStyle(1);;
ntuple->SetLineWidth(2); 
ntuple->Draw("nucZ", "typePart==0","same");

// Double differential (particle directions)
r->cd(5);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 

ntuple->SetLineColor(16);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("momZ", "typePart==7");
//ntuple->Draw("sqrt(momX^2+momY^2+momZ^2)", "typePart==7" ,"same");

ntuple->SetLineColor(19);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(2); 
ntuple->Draw("momZ", "typePart==3","same");
//ntuple->Draw("sqrt(momX^2+momY^2+momZ^2)", "typePart==7");

ntuple->SetLineColor(1);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(3); 
ntuple->Draw("momZ", "typePart==2", "same");
//ntuple->Draw("sqrt(momX^2+momY^2+momZ^2)", "typePart==1");
ntuple->SetLineColor(13);
ntuple->SetLineStyle(1);;
ntuple->SetLineWidth(2); 
ntuple->Draw("momZ", "typePart==2","same");
//ntuple->Draw("sqrt(momX^2+momY^2+momZ^2)", "typePart==2" ,"same");

r->cd(6);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 

ntuple->SetLineColor(19);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(4); 
//ntuple->Draw("sqrt(momX^2 + momY^2)/momZ +1", "typePart==7");
ntuple->SetLineColor(1);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("sqrt(momX*momX + momY*momY)", "typePart==1");

ntuple->SetLineColor(13);
ntuple->SetLineStyle(1);;
ntuple->SetLineWidth(2); 
//ntuple->Draw("momZ", "typePart==1");
ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==7","same");

ntuple->SetLineColor(16);

ntuple->SetLineStyle(1);
ntuple->SetLineWidth(3); 
//ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==3", "same");


r->Update();

r->Print("cascade.eps");

printf(" ::: Done.\n");
}



