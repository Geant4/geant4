{
// Batch analysis for g4 cascade simulation: root -b -q cascade.C
// Gray shades 10...19, 0 (from flack to white) 40-49 shades of red
// enum particleType { fragment = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

gROOT->Reset();

printf(" ::: Reading cascade data ...\n");

ifstream in;
in.open("cascade.out", ios::in);

Float_t nEve, typePart, eKin, momX, momY, momZ, nucEx;
Int_t nucA, nucZ;
Int_t nlines = 0;

TFile *f = new TFile("cascade.root","RECREATE");
TH1F *h1 = new TH1F("h1","particle types",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from cascade.out","nEve:typePart:eKin:momX:momY:momZ:nucA:nucZ:nucEx");

while (1) {
  in >> nEve >> typePart >> eKin >> momX >> momY >> momZ >> nucA >> nucZ >> nucEx;
  if (!in.good()) break;

  // Show first lines of data
  if (nlines < 2) printf("nEve = %1i, typePart = %1i, eKin = %4f\n",nEve, typePart, eKin);

  h1->Fill(typePart);
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
  gPad->SetGrid();
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
ntuple->SetFillColor(18);
ntuple->Draw("eKin", "typePart==1");

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

// pi+ 
ntuple->SetLineColor(16);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("eKin", "typePart==3");

// pi-
ntuple->SetLineStyle(1);
ntuple->SetLineColor(1);
ntuple->SetLineWidth(2); 
ntuple->Draw("eKin", "typePart==5", "same");

// pi0
//ntuple->SetLineStyle(4); // dots
ntuple->SetLineColor(13);
ntuple->SetLineWidth(1); 
ntuple->Draw("eKin", "typePart==7", "same");

// Fotons
r->cd(4);
ntuple->SetLineColor(1);
ntuple->SetFillColor(0);
ntuple->Draw("eKin", "typePart==0");

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
//ntuple->Draw("momZ", "typePart==2","same");
//ntuple->Draw("sqrt(momX^2+momY^2+momZ^2)", "typePart==2" ,"same");


r->cd(6);

ntuple->SetFillStyle(4000);
ntuple->SetFillColor(0); 

ntuple->SetLineColor(19);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(4); 
ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==7");
ntuple->SetLineColor(1);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(1); 
ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==1");

ntuple->SetLineColor(13);
ntuple->SetLineStyle(1);;
ntuple->SetLineWidth(2); 
//ntuple->Draw("momZ", "typePart==1","same");
ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==2", "same");

ntuple->SetLineColor(16);
ntuple->SetLineStyle(1);
ntuple->SetLineWidth(3); 
ntuple->Draw("sqrt(momX^2 + momY^2)/momZ", "typePart==3", "same");


r->Update();

printf(" ::: Writing ps files ...\n");
r->Print("cascade.eps");

printf(" ::: Done.\n");
}

