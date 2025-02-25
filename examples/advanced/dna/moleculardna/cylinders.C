{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(00000);

  system ("hadd -O -f molecular-dna.root molecular-dna_t*.root");

  c1 = new TCanvas("c1", "Damages", 120, 60, 1000, 1000);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillStyle(4000);
  c1->Divide(2,1);
  gPad->SetLeftMargin(0.13);

  Int_t SSBd, SSBi, SSBm, DSBd, DSBi, DSBm, DSBh;
  Int_t total_SSBd, total_SSBi, total_SSBm, total_DSBd, total_DSBi, total_DSBm, total_DSBh;
  Float_t total_SSBd2, total_SSBi2, total_SSBm2, total_DSBd2, total_DSBi2, total_DSBm2, total_DSBh2;
  Float_t mean_SSBd, mean_SSBi, mean_SSBm, mean_DSBd, mean_DSBi, mean_DSBm, mean_DSBh;

  Int_t SSB, SSBp, twoSSB, DSB, DSBp, DSBpp;
  Int_t total_SSB, total_SSBp, total_twoSSB, total_DSB, total_DSBp, total_DSBpp;
  Int_t total_SSB2, total_SSBp2, total_twoSSB2, total_DSB2, total_DSBp2,
    total_DSBpp2;
  Float_t mean_SSB, mean_SSBp, mean_twoSSB, mean_DSB, mean_DSBp, mean_DSBpp;

  total_SSBd = 0;
  total_SSBi = 0;
  total_SSBm = 0;
  total_DSBd = 0;
  total_DSBi = 0;
  total_DSBm = 0;
  total_DSBh = 0;

  total_SSBd2 = 0;
  total_SSBi2 = 0;
  total_SSBm2 = 0;
  total_DSBd2 = 0;
  total_DSBi2 = 0;
  total_DSBm2 = 0;
  total_DSBh2 = 0;

  total_SSB = 0;
  total_SSBp = 0;
  total_twoSSB = 0;
  total_DSB = 0;
  total_DSBp = 0;
  total_DSBpp = 0;

  total_SSB2 = 0;
  total_SSBp2 = 0;
  total_twoSSB2 = 0;
  total_DSB2 = 0;
  total_DSBp2 = 0;
  total_DSBpp2 = 0;

  Double_t EnergyDeposited_eV = 0;
  Double_t acc_edep = 0;
  Double_t acc_edep2 = 0;

  Double_t Energy;
  Char_t Primary;

  TFile *f = TFile::Open("molecular-dna.root");
  TTree *tree = (TTree *) f->Get("tuples/primary_source");
  Float_t number = (Float_t) tree->GetEntries();

  tree = (TTree *) f->Get("tuples/source");
  tree->SetBranchAddress("Primary",&Primary);
  tree->SetBranchAddress("Energy",&Energy);
  tree->SetBranchAddress("SSBd",&SSBd);
  tree->SetBranchAddress("SSBi",&SSBi);
  tree->SetBranchAddress("SSBm",&SSBm);
  tree->SetBranchAddress("DSBd",&DSBd);
  tree->SetBranchAddress("DSBi",&DSBi);
  tree->SetBranchAddress("DSBm",&DSBm);
  tree->SetBranchAddress("DSBh",&DSBh);

  Long64_t nentries = tree->GetEntries();

  for(int i = 0;i<nentries;i++){
    tree->GetEntry(i);
    total_SSBd += SSBd;
    total_SSBd2 += SSBd *SSBd;
    total_SSBi += SSBi;
    total_SSBi2 += SSBi *SSBi;
    total_SSBm += SSBm;
    total_SSBm2 += SSBm *SSBm;

    total_DSBd += DSBd;
    total_DSBd2 += DSBd *DSBd;
    total_DSBi += DSBi;
    total_DSBi2 += DSBi *DSBi;
    total_DSBm += DSBm;
    total_DSBm2 += DSBm *DSBm;
    total_DSBh += DSBh;
    total_DSBh2 += DSBh *DSBh;
  }

  tree = (TTree *) f->Get("tuples/damage");
  tree->SetBranchAddress("EnergyDeposited_eV",&EnergyDeposited_eV);
  nentries = tree->GetEntries();
  for(int i = 0;i<nentries;i++){
    tree->GetEntry(i);
    acc_edep += EnergyDeposited_eV;
    acc_edep2 += EnergyDeposited_eV *EnergyDeposited_eV;
  }

  tree = (TTree *) f->Get("tuples/classification");
  tree->SetBranchAddress("SSB",&SSB);
  tree->SetBranchAddress("SSBp",&SSBp);
  tree->SetBranchAddress("2SSB",&twoSSB);
  tree->SetBranchAddress("DSB",&DSB);
  tree->SetBranchAddress("DSBp",&DSBp);
  tree->SetBranchAddress("DSBpp",&DSBpp);

  nentries = tree->GetEntries();

  for(int i = 0;i<nentries;i++){
    tree->GetEntry(i);
    total_SSB += SSB;
    total_SSB2 += SSB *SSB;
    total_SSBp += SSBp;
    total_SSBp2 += SSBp *SSBp;
    total_twoSSB += twoSSB;
    total_twoSSB2 += twoSSB *twoSSB;

    total_DSB += DSB;
    total_DSB2 += DSB *DSB;
    total_DSBp += DSBp;
    total_DSBp2 += DSBp *DSBp;
    total_DSBpp += DSBpp;
    total_DSBpp2 += DSBpp *DSBpp;
  }

  mean_SSBd = (Float_t) total_SSBd / number;
  mean_SSBi = (Float_t) total_SSBi / number;
  mean_SSBm = (Float_t) total_SSBm / number;

  mean_DSBd = (Float_t) total_DSBd / number;
  mean_DSBi = (Float_t) total_DSBi / number;
  mean_DSBm = (Float_t) total_DSBm / number;
  mean_DSBh = (Float_t) total_DSBh / number;

  Double_t SD_SSBd = sqrt(((total_SSBd2 / number) - pow(total_SSBd / number,2))/(number -1));
  Double_t SD_SSBi = sqrt(((total_SSBi2 / number) - pow(total_SSBi / number,2))/(number -1));
  Double_t SD_SSBm = sqrt(((total_SSBm2 / number) - pow(total_SSBm / number,2))/(number -1));

  Double_t SD_DSBd = sqrt(((total_DSBd2 / number) - pow(total_DSBd / number,2))/(number -1));
  Double_t SD_DSBi = sqrt(((total_DSBi2 / number) - pow(total_DSBi / number,2))/(number -1));
  Double_t SD_DSBm = sqrt(((total_DSBm2 / number) - pow(total_DSBm / number,2))/(number -1));
  Double_t SD_DSBh = sqrt(((total_DSBh2 / number) - pow(total_DSBh / number,2))/(number -1));


  mean_SSB = (Float_t) total_SSB / number;
  mean_SSBp = (Float_t) total_SSBp / number;
  mean_twoSSB = (Float_t) total_twoSSB / number;

  mean_DSB = (Float_t) total_DSB / number;
  mean_DSBp = (Float_t) total_DSBp / number;
  mean_DSBpp = (Float_t) total_DSBpp / number;

  Double_t SD_SSB = sqrt(((total_SSB2 / number) - pow(total_SSB / number,2))/(number -1));
  Double_t SD_SSBp = sqrt(((total_SSBp2 / number) - pow(total_SSBp / number,2))
                            /(number -1));
  Double_t SD_twoSSB = sqrt(((total_twoSSB2 / number) - pow(total_twoSSB /
                                                              number,2))
                            /(number -1));

  Double_t SD_DSB = sqrt(((total_DSB2 / number) - pow(total_DSB / number,2))/(number -1));
  Double_t SD_DSBp = sqrt(((total_DSBp2 / number) - pow(total_DSBp / number,2))
                           /(number -1));
  Double_t SD_DSBpp = sqrt(((total_DSBpp2 / number) - pow(total_DSBpp / number,2))
                           /(number -1));

  cout<<"Paricle : "<<Primary<<'\t'
       <<"Energy [/MeV] : "<<Energy<<'\t'
       <<"number : "<<number<<'\n'
       <<" Output Damages : "<<'\n'<<
    '\t'<<"SSB direct : "<<mean_SSBd<<'\t'<<" error : "<<SD_SSBd<<'\n'<<
    '\t'<<"SSB indirect : "<<mean_SSBi<<'\t'<<" error : "<<SD_SSBi<<'\n'<<
    '\t'<<"SSB mix : "<<mean_SSBm<<'\t'<<" error : "<<SD_SSBm<<'\n'<<'\n'<<

    '\t'<<"DSB direct : "<<mean_DSBd<<'\t'<<" error : "<<SD_DSBd<<'\n'<<
    '\t'<<"DSB indirect : "<<mean_DSBi<<'\t'<<" error : "<<SD_DSBi<<'\n'<<
    '\t'<<"DSB mix : "<<mean_DSBm<<'\t'<<" error : "<<SD_DSBm<<'\n'<<
    '\t'<<"DSB hybrid : "<<mean_DSBh<<'\t'<<" error : "<<SD_DSBh<<'\n'<<'\n'<<

    '\t'<<"SSB  : "<<mean_SSB<<'\t'<<" error : "<<SD_SSB<<'\n'<<
    '\t'<<"SSBp : "<<mean_SSBp<<'\t'<<" error : "<<SD_SSBp<<'\n'<<
    '\t'<<"2SSB : "<<mean_twoSSB<<'\t'<<" error : "<<SD_twoSSB<<'\n'<<'\n'<<

    '\t'<<"DSB : "<<mean_DSB<<'\t'<<" error : "<<SD_DSB<<'\n'<<
    '\t'<<"DSBp : "<<mean_DSBp<<'\t'<<" error : "<<SD_DSBp<<'\n'<<
    '\t'<<"DSBpp : "<<mean_DSBpp<<'\t'<<" error : "<<SD_DSBpp<<'\n';


  f->Close();
  const Int_t n1 = 7;
  Double_t x1[n1] = {2,4,6,8,10,12,14};
  Double_t _y1[n1] = {mean_SSBd,mean_SSBi,mean_SSBm,mean_DSBd,mean_DSBi,
                       mean_DSBm,mean_DSBh};
  Double_t err_y1[n1] = {SD_SSBd,SD_SSBi,SD_SSBm,SD_DSBd,SD_DSBi,SD_DSBm,SD_DSBh};
  TGraph* gr1 = new TGraphErrors(n1,x1,_y1,0,err_y1);
  gr1->SetTitle("Break Source Frequency");
  gr1->GetXaxis()->SetBinLabel(10,"SSB direct");
  gr1->GetXaxis()->SetBinLabel(20,"SSB indirect");
  gr1->GetXaxis()->SetBinLabel(35,"SSB mix");
  gr1->GetXaxis()->SetBinLabel(50,"DSB direct");
  gr1->GetXaxis()->SetBinLabel(65,"DSB indirect");
  gr1->GetXaxis()->SetBinLabel(80,"DSB mix");
  gr1->GetXaxis()->SetBinLabel(90,"DSB hybrid");
  gr1->SetFillColor(49);

  const Int_t n2 = 6;
  Double_t x2[n2] = {2,4,6,8,10,12};
  Double_t _y2[n2] = {mean_SSB,mean_SSBp,mean_twoSSB,mean_DSB,mean_DSBp,
                       mean_DSBpp};
  Double_t err_y2[n2] = {SD_SSB,SD_SSBp,SD_twoSSB,SD_DSB,SD_DSBp,SD_DSBpp};
  TGraph* gr2 = new TGraphErrors(n2,x2,_y2,0,err_y2);
  gr2->SetTitle("Break Complexity Frequency");
  gr2->GetXaxis()->SetBinLabel(10,"SSB");
  gr2->GetXaxis()->SetBinLabel(25,"SSBp");
  gr2->GetXaxis()->SetBinLabel(45,"2SSB");
  gr2->GetXaxis()->SetBinLabel(60,"DSB");
  gr2->GetXaxis()->SetBinLabel(75,"DSBp");
  gr2->GetXaxis()->SetBinLabel(90,"DSBpp");
  gr2->SetFillColor(49);

  //Draw
  c1->cd(1);
  gr1->Draw("ba");

  c1->cd(2);
  gr2->Draw("ba");
}
