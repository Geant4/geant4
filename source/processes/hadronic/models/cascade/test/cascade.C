{
// Batch analysis for g4 cascade simulation: root -b -q cascade.C

#include <iostream.h>


gROOT->Reset();

// enum particleType { proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

ifstream in;
in.open("cascade.out", ios::in);

Float_t nEve, typePart, eKin, momX, momY, momZ;
Int_t nlines = 0;

TFile *f = new TFile("cascade.root","RECREATE");
TH1F *h1 = new TH1F("h1","particle types",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from cascade.out","nEve:typePart:eKin:momX:momY:momZ");

printf(" ::: Reading cascade data ...\n");
while (1) {
  in >> nEve >> typePart >> eKin >> momX >> momY >> momZ;
  if (!in.good()) break;
  if (nlines < 3) printf("nEve = %1i, typePart = %1i, eKin = %4f\n",nEve, typePart, eKin);
  h1->Fill(typePart);
  ntuple->Fill(nEve, typePart, eKin, momX, momY, momZ);
  nlines++;
}

printf(" ::: Found %d points\n",nlines);

in.close();

f->Write();

c1 = new TCanvas("c1","The Ntuple canvas",200,10,700,780);

pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.48,0.98,21);
pad2 = new TPad("pad2","This is pad2",0.52,0.52,0.98,0.98,21);
pad3 = new TPad("pad3","This is pad3",0.02,0.02,0.48,0.48,21);
pad4 = new TPad("pad4","This is pad4",0.52,0.02,0.98,0.48,1);
pad1->Draw();
pad2->Draw();
pad3->Draw();
pad4->Draw();

// Change default style for the statistics box
gStyle->SetStatW(0.30);
gStyle->SetStatH(0.20);
gStyle->SetStatColor(42);


// Proton energy
pad1->cd();
pad1->SetGrid();
pad1->SetLogy();
pad1->GetFrame()->SetFillColor(15);
ntuple->SetLineColor(1);
ntuple->SetFillStyle(1001);
ntuple->SetFillColor(45);

ntuple->Draw("eKin", "typePart==1");
c1->Update();

// Neutron energy
pad2->cd();
pad2->SetGrid();
pad2->SetLogy();
pad2->GetFrame()->SetFillColor(15);
ntuple->SetLineColor(1);
ntuple->SetFillStyle(1001);
ntuple->SetFillColor(45);

ntuple->Draw("eKin", "typePart==2");
c1->Update();


// Display a scatter plot of two columns with a selection.
// Superimpose the result of another cut with a different marker color
pad3->cd();
pad3->SetGrid();
pad3->SetLogy();
pad3->GetFrame()->SetFillColor(15);
ntuple->SetLineColor(1);
ntuple->SetFillStyle(3001);
ntuple->SetFillColor(45);

ntuple->Draw("eKin", "typePart==3");
ntuple->Draw("eKin", "typePart==5");
ntuple->Draw("eKin", "typePart==7");
c1->Update();

// Display a 3-D scatter plot of 3 columns. Superimpose a different selection.
pad4->cd();
pad4->SetGrid();
pad4->SetLogy();
pad4->GetFrame()->SetFillColor(15);
ntuple->SetLineColor(1);
ntuple->SetFillStyle(1001);
ntuple->SetFillColor(45);

ntuple->Draw("eKin", "typePart==10");
c1->Update();

c1->cd();
c1->Update();
gStyle->SetStatColor(19);

printf(" ::: Writing ps files ...\n");
c1->Print("cascade.ps");

printf(" ::: Done.\n");
}


