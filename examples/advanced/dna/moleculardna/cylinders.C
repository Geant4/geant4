{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(00000);

  auto c1 = new TCanvas("c1", "Damages", 120, 60, 1000, 1000);
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
  //Double_t acc_edep = 0;
  //Double_t acc_edep2 = 0;

  Double_t Energy;
  Char_t Primary;

  TFile *f = TFile::Open("molecular-dna.root");
  TTree *tree = (TTree *) f->Get("tuples/primary_source");
  Long64_t number = (Float_t) tree->GetEntries();

  if (number<2) {
    std::cout << "Not enough entries in the \"primary_source\" TTree (" << (long)number << " entries)\n";
    gApplication->Terminate(0);
  }

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

  Long64_t nentriesS = tree->GetEntries();

  for(int i = 0;i<nentriesS;i++){
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

  /*
  tree = (TTree *) f->Get("tuples/damage");
  tree->SetBranchAddress("EnergyDeposited_eV",&EnergyDeposited_eV);
  nentries = tree->GetEntries();
  for(int i = 0;i<nentries;i++){
    tree->GetEntry(i);
    acc_edep += EnergyDeposited_eV;
    acc_edep2 += EnergyDeposited_eV *EnergyDeposited_eV;
  }
  */

  tree = (TTree *) f->Get("tuples/classification");
  tree->SetBranchAddress("SSB",&SSB);
  tree->SetBranchAddress("SSBp",&SSBp);
  tree->SetBranchAddress("2SSB",&twoSSB);
  tree->SetBranchAddress("DSB",&DSB);
  tree->SetBranchAddress("DSBp",&DSBp);
  tree->SetBranchAddress("DSBpp",&DSBpp);

  Long64_t nentriesC = tree->GetEntries();

  for(int i = 0;i<nentriesC;i++){
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

  // Mean values 

  mean_SSBd = (Float_t) total_SSBd / nentriesS;
  mean_SSBi = (Float_t) total_SSBi / nentriesS;
  mean_SSBm = (Float_t) total_SSBm / nentriesS;

  mean_DSBd = (Float_t) total_DSBd / nentriesS;
  mean_DSBi = (Float_t) total_DSBi / nentriesS;
  mean_DSBm = (Float_t) total_DSBm / nentriesS;
  mean_DSBh = (Float_t) total_DSBh / nentriesS;

  // SEM values

  Double_t SD_SSBd = sqrt(abs(((total_SSBd2 / nentriesS) - pow(total_SSBd / nentriesS,2)))/(nentriesS -1));
  Double_t SD_SSBi = sqrt(abs(((total_SSBi2 / nentriesS) - pow(total_SSBi / nentriesS,2)))/(nentriesS -1));
  Double_t SD_SSBm = sqrt(abs(((total_SSBm2 / nentriesS) - pow(total_SSBm / nentriesS,2)))/(nentriesS -1));

  Double_t SD_DSBd = sqrt(abs(((total_DSBd2 / nentriesS) - pow(total_DSBd / nentriesS,2)))/(nentriesS -1));
  Double_t SD_DSBi = sqrt(abs(((total_DSBi2 / nentriesS) - pow(total_DSBi / nentriesS,2)))/(nentriesS -1));
  Double_t SD_DSBm = sqrt(abs(((total_DSBm2 / nentriesS) - pow(total_DSBm / nentriesS,2)))/(nentriesS -1));
  Double_t SD_DSBh = sqrt(abs(((total_DSBh2 / nentriesS) - pow(total_DSBh / nentriesS,2)))/(nentriesS -1));

  // Mean values
  
  mean_SSB = (Float_t) total_SSB / nentriesC;
  mean_SSBp = (Float_t) total_SSBp / nentriesC;
  mean_twoSSB = (Float_t) total_twoSSB / nentriesC;

  mean_DSB = (Float_t) total_DSB / nentriesC;
  mean_DSBp = (Float_t) total_DSBp / nentriesC;
  mean_DSBpp = (Float_t) total_DSBpp / nentriesC;

  // SEM values
  
  Double_t SD_SSB = sqrt(abs(((total_SSB2 / nentriesC) - pow(total_SSB / nentriesC,2))/(nentriesC -1)));
  Double_t SD_SSBp = sqrt(abs(((total_SSBp2 / nentriesC) - pow(total_SSBp / nentriesC,2))
                            /(nentriesC -1)));
  Double_t SD_twoSSB = sqrt(abs(((total_twoSSB2 / nentriesC) - pow(total_twoSSB /
                                                              nentriesC,2))
                            /(nentriesC -1)));

  Double_t SD_DSB = sqrt(abs(((total_DSB2 / nentriesC) - pow(total_DSB / nentriesC,2))/(nentriesC -1)));
  Double_t SD_DSBp = sqrt(abs(((total_DSBp2 / nentriesC) - pow(total_DSBp / nentriesC,2))
                           /(nentriesC -1)));
  Double_t SD_DSBpp = sqrt(abs(((total_DSBpp2 / nentriesC) - pow(total_DSBpp / nentriesC,2))
                           /(nentriesC -1)));

  //
  
  cout<<"Particle : "<<Primary<<'\t'
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

  c1->cd(1);
  gr1->Draw("ba");

  c1->cd(2);
  gr2->Draw("ba");
}
