Double_t expdist(Double_t *x, Double_t *par)
{
  Double_t fitval = par[0]*TMath::Exp(-1.0*x[0]);
   return fitval;
}
void myfit()
{
  gROOT->Reset();
  #include <iostream.h>
  Float_t x;

  //gBenchmark->Start("fit1");
  TH1F *h1 = new TH1F("h1","exprnf",40,0,6);

  //cout << "Here we are 0" << endl;
  
  ifstream in;
  in.open("exprnf.out", ios::in);
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
   TF1 *func = new TF1("expdist",expdist,0,6,1);
   
   // cout << "Here we are 3" << endl; 

   // Sets initial values and parameter names
   // func->SetParameters(1.0);
   //func->SetParNames("Constant");

   //cout << "Here we are 4" << endl; 
// Fit histogram in range defined by function
   h1->Fit("expdist", "r");
   h1->GetXaxis()->SetTitle("return value");
   h1->GetXaxis()->CenterTitle();
   h1->GetYaxis()->SetTitle("number of hits");
   h1->GetYaxis()->CenterTitle();

   //gBenchmark->Show("fit1");
// Gets integral of function between fit limits
//   printf("Integral of function = %g\n",func->Integral(-2,2));
}


