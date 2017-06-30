// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
	
  c1 = new TCanvas ("c1","",20,20,800,800);
  c1->Divide(1,1);

  TFile f("slowing.root"); 

  TH1F* h1 ;
  h1 = (TH1F*)f.Get("1"); 
  h2 = (TH1F*)f.Get("2"); 
  h3 = (TH1F*)f.Get("3"); 
     
//goto end;

  Int_t nbinsx = h1->GetXaxis()->GetNbins();
  //cout << nbinsx << endl;

  Double_t y = 0;
  Double_t mini = 0;
  Double_t maxi = 0;
  Double_t largeur = 0;

  Double_t sum = 0;

  // Division by bin width to get y axis 
  // in nm/eV
  //
  // Scaling by 1E9/1.6 to get correct unit 
  // for Phi/D in (/cm2/eV/Gy) 
  // when histogram (in nm/eV) is
  // multiplied by density(=1g/cm3)/E(eV)

  for (Int_t i=1; i<=nbinsx; i++)
  {
    sum = sum + h1->GetBinContent(i);

    mini = h1->GetBinLowEdge(i);
    maxi = mini + h1->GetBinWidth(i);
    largeur = std::pow(10,maxi)-std::pow(10,mini);
    // cout << mini << " " << std::pow(10,mini)<< " " << largeur 
    // << " " << maxi << " " << std::pow(10,maxi) << endl;
    h1->SetBinContent(i,h1->GetBinContent(i)*(1E9/1.6)/largeur);
    h2->SetBinContent(i,h2->GetBinContent(i)*(1E9/1.6)/largeur);
    h3->SetBinContent(i,h3->GetBinContent(i)*(1E9/1.6)/largeur);

  }

  gStyle->SetOptStat(000000);

  cout << endl;
  cout << "--> Integral of Phi (nm/eV) = " << sum << endl;
  cout << endl;

c1->cd(1);

  TH2F *ht = new TH2F("","",2,1,6,2,1E2,1E8);
  ht->Draw();
  ht->GetXaxis()->SetTitle("Log(E (eV))");
  ht->GetYaxis()->SetTitle("#phi/D (/cm^{2}/eV/Gy)");
  ht->GetXaxis()->SetTitleSize(0.03);
  ht->GetYaxis()->SetTitleSize(0.03);
  ht->GetXaxis()->SetTitleOffset(1.7);
  ht->GetYaxis()->SetTitleOffset(1.7);

  gPad->SetLogy();
  h1->SetLineColor(2);
  h1->Draw("HSAME");
  h2->SetLineColor(3);
  h2->Draw("HSAME");
  h3->SetLineColor(4);
  h3->Draw("HSAME");
  h1->Draw("HSAME");

  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->AddEntry(h1,"All e-","L");
  legend->AddEntry(h2,"Primaries","L");
  legend->AddEntry(h3,"Secondaries","L");
  legend->Draw();

end:
}
