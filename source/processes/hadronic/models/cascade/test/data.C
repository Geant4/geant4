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

//printf(" ::: Writing ps file ...\n");
c1->Print("data.eps");

printf(" ::: Done.\n");
}

