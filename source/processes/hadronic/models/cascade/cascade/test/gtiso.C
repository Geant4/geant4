Double_t unidist(Double_t *x, Double_t *par)
{
  Double_t fitval = par[0];
   return fitval;
}

void myfit()
{
  gROOT->Reset();
  #include <iostream.h>
  Float_t x;

  //gBenchmark->Start("fit1");
  TH1F *h1 = new TH1F("h1","gtiso",60,-1,1);

  //cout << "Here we are 0" << endl;
  
  ifstream in;
  in.open("gtiso.out", ios::in);
  while (1) {
  in >> x;
  if (!in.good()) break;
  h1->Fill(x);
 }  
  in.close();

  //cout << "Here we are" << endl;

  TCanvas *c1 = new TCanvas("c1","the fit canvas",500,400);
  c1->SetFillColor(0);

  //cout << "Here we are 2" << endl;

// Creates a Root function based on function fitf above
   TF1 *func = new TF1("unidist",unidist,-1,1,1);
   
   // cout << "Here we are 3" << endl; 

// Fit histogram in range defined by function
   h1->Fit("unidist", "r");
   h1->GetXaxis()->SetTitle("parameter 'u'");
   h1->GetXaxis()->CenterTitle();
   h1->GetYaxis()->SetTitle("number of hits");
   h1->GetYaxis()->CenterTitle();

   //gBenchmark->Show("fit1");
// Gets integral of function between fit limits
//   printf("Integral of function = %g\n",func->Integral(-2,2));
}


