{
// Batch analysis for g4 cascade data : root -b -q data.C

gROOT->Reset();

printf(" ::: Reading cascade data ...\n");

ifstream in;
in.open("data.out", ios::in);

Float_t e[100], xc[8][100];
Int_t nlines = 0;

for (Int_t j = 1; j < 100; j++) {

  in >> e[j];
  printf("j = %1i, e = %4f\n",j, e[j]);

  for (Int_t i = 1; i < 7; i++) {
    in >> xc[i][j];

    printf("i = %1i, xc = %4f\n",i, xc[i][j]);
  };
}

c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
c1->SetFillColor(0);
c1->SetBorderMode(0);
c1->Divide(2,2);
  c1->cd(1);
  gPad->SetBorderMode(0);
c1->SetFillColor(11);
c1->SetGrid();

const Int_t n = 100;
Double_t x[n], y[n];

for (Int_t i=0;i<n;i++) {
  x[i] = e[i];
  y[i] = xc[1][i];
  printf(" e =  %i %f %f \n", i, x[i], y[i]);
}

gr = new TGraph(n,x,y);
gr->SetLineColor(1);
gr->SetLineWidth(1);
gr->SetMarkerColor(17);
gr->SetMarkerStyle(21);
gr->SetTitle("Cascade cross section data");
gr->Draw("ACP");

c1->Update();
c1->GetFrame()->SetFillColor(13);
c1->GetFrame()->SetBorderSize(0);
gr->GetHistogram()->SetXTitle("GeV");
gr->GetHistogram()->SetYTitle("mb");
c1->Modified();

// ------------------------------------------
Double_t t[5][5] = 
{{0.5, 0.51, 2.37, 2.97, 4.83 },
 {1.0, 0.38, 2.70, 3.87, 6.95 },
 {2.0, 0.54, 3.02, 5.12, 10.50},
 {4.0, 0.46, 3.78, 6.07, 14.91},
 {8.0, 0.43, 4.12, 7.06, 18.95}};

c1->cd(2);
  gPad->SetBorderMode(0);
// gPad->SetLogy();
c1->SetFillColor(0);
c1->SetGrid();

const Int_t n = 5;
Double_t x[n], tH[n], tAl[n], tFe[n], tPb[n];

for (Int_t i=0; i < n; i++) {
    x[i] = t[i][0];
   tH[i] = t[i][1];
  tAl[i] = t[i][2];
  tFe[i] = t[i][3];
  tPb[i] = t[i][4];
  printf(" e =  %i %f %f \n", i, x[i], tH[i]);
}


gr0 = new TGraph(n, x, tPb);
gr1 = new TGraph(n, x, tFe);
gr2 = new TGraph(n, x, tAl);
gr3 = new TGraph(n, x, tH);


gr0->SetLineColor(1);
gr0->SetLineWidth(1);

gr1->SetLineColor(1);
gr1->SetLineWidth(1);

gr2->SetLineColor(1);
gr2->SetLineWidth(1);

gr3->SetLineColor(1);
gr3->SetLineWidth(1);

gr0->SetMarkerColor(16);
gr1->SetMarkerColor(16);
gr2->SetMarkerColor(16);
gr3->SetMarkerColor(16);

gr0->SetMarkerStyle(21);
gr1->SetMarkerStyle(21);
gr2->SetMarkerStyle(21);
gr3->SetMarkerStyle(21);

gr0->SetTitle("Simulation time of one INC");

gr0->Draw("AC");
gr1->Draw("C");
gr2->Draw("C");
gr3->Draw("C");
TLatex l;
l->DrawLatex(0.5,  19, "Geant4 INC model speed");
l->DrawLatex(5,  18, "Pb_{82}^{208}");
l->DrawLatex(5, 7.5, "Fe_{26}^{56}");
l->DrawLatex(5, 4.5, "Al_{13}^{26}");

c1->Update();
c1->GetFrame()->SetFillColor(0);
c1->GetFrame()->SetBorderSize(0);

gr0->GetHistogram()->SetXTitle("Bullet particle kinetic energy [GeV]");
gr0->GetHistogram()->SetYTitle("Timing for one INC [ms] ");

c1->Modified();

c1->Print("data.eps");

printf(" ::: Done.\n");
}

