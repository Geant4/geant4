//
// Fit of longitudinal profile with Gamma function
// See Review of Particles Physics - Electromagnetic cascades
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Double_t GammaFunction(Double_t* t, Double_t* params ) {
  if (ROOT::TMath::Gamma(params[0]) != 0 ) {
   return  (100 * pow(params[1],params[0]) / ROOT::TMath::Gamma(params[0]))
           * pow(t[0],(params[0]-1)) * exp (-params[1]*t[0]);
  } else {	 
   return 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GammaFit() {
 cout 
  << "Try fit [100*pow(b,a)/ROOT::TMath::Gamma(a)] * pow(t,(a-1)) * exp(-b*t)"  
  << endl;
TFile* f = new  TFile("93ref0.root");
TH1D* h1 = (TH1D*) f->Get("4");
//put error in h1
TH1D* h2 = (TH1D*) f->Get("5");
int nmax = h1->GetNbinsX();
for (int n=0; n<nmax; n++) {
  Double_t er = h2->GetBinContent(n);
  h1->SetBinError(n,er);
}

TF1* func = new TF1("fit", GammaFunction,-0.5,40.,2);
func->SetParameter(0,1.);
func->SetParameter(1,1.);
func->SetParNames("a", "b");

h1->Fit("fit","r","");
h1->SetMarkerColor(kRed);
h1->SetMarkerStyle(3);
h1->Draw("epl");

Double_t a = func->GetParameter(0);
Double_t b = func->GetParameter(1);
Double_t tmax = (a-1)/b;

// Chi-square distribution
double chisq=func->GetChisquare(); 
double ndf=func->GetNDF();
double chisqdf=chisq/ndf;

gStyle->SetOptFit(1111);

cout << "----------------------------------------------";
cout << "\n Chisquare: " << chisq << " / " << ndf << " = " << chisqdf;
cout << "\n----------------------------------------------";
cout << "\n a = " << a;
cout << "\n b = " << b;
cout << "\n tmax = " << tmax;
cout << "\n----------------------------------------------" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
