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

c1 = new TCanvas("c1","A Simple Graph Example", 200, 10, 700, 500);
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
gr->SetMarkerStyle(1);
gr->SetTitle("");
gr->Draw("ACP");

c1->Update();
c1->GetFrame()->SetFillColor(13);
c1->GetFrame()->SetBorderSize(0);
gr->GetHistogram()->SetXTitle("Particle energy [GeV]");
gr->GetHistogram()->SetYTitle("cross-section mb");
TLatex l;
l->DrawLatex(0.5,  65, "Cross-section p-:::");
l->DrawLatex(5,  18, "Pb");
c1->Modified();

// Timing plot
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

gr0->SetTitle("");

gr0->Draw("AC");
gr1->Draw("C");
gr2->Draw("C");
gr3->Draw("C");
TLatex l;
l->DrawLatex(0.5,  19, "Geant4 INC model speed");
l->DrawLatex(5,  18, "Pb");
l->DrawLatex(5, 7.5, "Fe");
l->DrawLatex(5, 4.5, "Al");

c1->Update();
c1->GetFrame()->SetFillColor(0);
c1->GetFrame()->SetBorderSize(0);

gr0->GetHistogram()->SetXTitle("Kinetic energy of bullet [GeV]");
gr0->GetHistogram()->SetYTitle("Timing for one INC [ms] ");

c1->Modified();


// particle production crossection
Double_t cs[31][6] = 
{{   5, 1.3906, 0.6094, 1.0888, 0.4794,        0.776}, 
 { 4.5, 1.379 , 0.621 , 1.0973, 0.4763,        0.734}, 
 {   4, 1.3044, 0.6956, 0.9944, 0.2988,       0.5662}, 
 { 3.5, 1.2962, 0.7038, 0.9334, 0.2296,        0.536}, 
 {   3, 1.2956, 0.7044, 0.9068, 0.2024,       0.5134}, 
 { 2.5, 1.2928, 0.7072, 0.8494, 0.1422,       0.4641}, 
 {   2, 1.2848, 0.7152, 0.7975, 0.0823,       0.3998}, 
 { 1.9, 1.2754, 0.7246, 0.7985, 0.0739,       0.3734}, 
 { 1.8, 1.2751, 0.7249, 0.7921, 0.0672,        0.354},  
 { 1.7, 1.2703, 0.7297, 0.7869, 0.0572,        0.335},  
 { 1.6, 1.2586, 0.7414, 0.7868, 0.0454,       0.3178},  
 { 1.5, 1.2471, 0.7529, 0.7826, 0.0297,        0.298},  
 { 1.4, 1.2352, 0.7648, 0.7887, 0.0239,       0.2737},  
 { 1.3, 1.2282, 0.7718, 0.7859, 0.0141,       0.2619},  
 { 1.2, 1.2169, 0.7831, 0.7893, 0.0062,       0.2355},  
 { 1.1, 1.2137, 0.7863, 0.7883, 0.002 ,       0.2244},  
 {   1, 1.2051, 0.7949, 0.7955, 0.0006,       0.2084},  
 { 0.9, 1.1939, 0.8061, 0.8062, 0.0001,        0.195},  
 {0.85, 1.1904, 0.8096, 0.8096,      0,       0.1908},  
 { 0.8, 1.1863, 0.8137, 0.8137,      0,       0.1863},  
 {0.75, 1.1809, 0.8191, 0.8191,      0,        0.181},  
 { 0.7, 1.1756, 0.8244, 0.8244,      0,       0.1756},  
 {0.65, 1.1719, 0.8281, 0.8281,      0,       0.1719},  
 { 0.6, 1.1563, 0.8437, 0.8437,      0,       0.1563},  
 {0.585,1.1533, 0.8467, 0.8467,      0,       0.1533},  
 {0.55, 1.1482, 0.8518, 0.8518,      0,       0.1482},  
 { 0.5, 1.131 , 0.869 , 0.869 ,      0,        0.131},  
 {0.45, 1.1332, 0.8668, 0.8668,      0,       0.1332},  
 { 0.4, 1.1273, 0.8727, 0.8727,      0,       0.1273},  
 {0.35, 1.1529, 0.8471, 0.8471,      0,       0.1529}};  


c1->cd(3);
gPad->SetBorderMode(0);
c1->SetFillColor(0);
c1->SetGrid();

const Int_t n = 31;
Double_t eKin[n], np[n], nn[n], npip[n], npim[n], npiz[n], nTot[n];

for (Int_t i = 0; i < n; i++) {
    eKin[i] = cs[i][0];
   np[i] = cs[i][1];
  nn[i] = cs[i][2];
  npip[i] = cs[i][3];
  npim[i] = cs[i][4];
  npiz[i] = cs[i][5];
  nTot[i] = cs[i][1] + cs[i][2] + cs[i][3] + cs[i][4] + cs[i][5];
  printf(" line =  %i  energy [GeV] = %f number of pi- = %f \n", i, eKin[i], npim[i]);
}
gr0 = new TGraph(n, eKin, np);
gr1 = new TGraph(n, eKin, nn);
gr2 = new TGraph(n, eKin, npip);
gr3 = new TGraph(n, eKin, npim);
gr4 = new TGraph(n, eKin, npiz);

gr0->SetLineColor(1);
gr0->SetLineWidth(1);

gr1->SetLineColor(1);
gr1->SetLineWidth(1);

gr2->SetLineColor(1);
gr2->SetLineWidth(1);

gr3->SetLineColor(1);
gr3->SetLineWidth(1);

gr4->SetLineColor(1);
gr4->SetLineWidth(1);


gr0->SetMarkerColor(16);
gr1->SetMarkerColor(16);
gr2->SetMarkerColor(16);
gr3->SetMarkerColor(16);
gr4->SetMarkerColor(16);

gr0->SetMarkerStyle(21);
gr1->SetMarkerStyle(21);
gr2->SetMarkerStyle(21);
gr3->SetMarkerStyle(21);
gr4->SetMarkerStyle(21);

gr0->SetTitle("");

gr0->Draw("AC");
gr1->Draw("C");
gr2->Draw("C");
gr3->Draw("C");
gr4->Draw("C");
TLatex l;
l->DrawLatex(0.1,  1.4, "10k collisions (p, H_{1}^{1})");
l->DrawLatex(5, 1.3, "p      ");
l->DrawLatex(5, 1.0, "pi ^{+}");
l->DrawLatex(5, 0.8, "pi ^{0}");
l->DrawLatex(5, 0.6, "n      ");
l->DrawLatex(5, 0.5, "pi ^{-}");

c1->Update();
c1->GetFrame()->SetFillColor(0);
c1->GetFrame()->SetBorderSize(0);

gr0->GetHistogram()->SetXTitle("Kinetic energy of bullet [GeV]");
gr0->GetHistogram()->SetYTitle("Number of particles");


c1->cd(4);
gPad->SetBorderMode(0);
c1->SetFillColor(0);
c1->SetGrid();

gr0 = new TGraph(n, eKin, nTot);

gr0->SetLineColor(1);
gr0->SetLineWidth(1);

gr0->SetMarkerColor(16);

gr0->SetMarkerStyle(21);

gr0->SetTitle("");

gr0->Draw("AC");

TLatex l;
l->DrawLatex(0.1,  4.4, "Total particle flux 10k (p, H_{1}^{1})");
//l->DrawLatex(5, 1.3, "p      ");


c1->Update();
c1->GetFrame()->SetFillColor(0);
c1->GetFrame()->SetBorderSize(0);

gr0->GetHistogram()->SetXTitle("Kinetic energy of bullet [GeV]");
gr0->GetHistogram()->SetYTitle("Number of particles");


c1->Modified();

c1->Print("data.eps");

printf(" ::: Done.\n");
}






