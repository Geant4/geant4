// *********************************************************************
// To execute this macro under ROOT after your simulation ended,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plot()
{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");

  TCanvas* c1 = new TCanvas ("c1","",20,20,1500,500);
  c1->Divide(3,1);

  // Uncomment if merging should be done
  //system ("rm -rf psp.root");
  //system ("hadd psp.root psp_*.root");

  TFile* f = new TFile("psp.root");

  TNtuple* ntuple;
  ntuple = (TNtuple*)f->Get("scorer");
  bool rowWise = true;
  TBranch* eventBranch = ntuple->FindBranch("row_wise_branch");
  if ( ! eventBranch ) rowWise = false;
  // std::cout <<  "rowWise: " << rowWise << std::endl;

  //*********************************************************************
  // canvas tab 1
  //*********************************************************************

  c1->cd(1);
  gStyle->SetOptStat(10);
  ntuple->SetFillStyle(1001);
  ntuple->SetFillColor(2);
  ntuple->Draw("PDGCode","","B");
  gPad->SetLogy();

  //*********************************************************************
  // canvas tab 2
  //*********************************************************************

  c1->cd(2);
  ntuple->SetMarkerColor(4);
  ntuple->SetMarkerSize(5);
  ntuple->Draw("x:y:z");

  //*********************************************************************
  // canvas tab 3
  //*********************************************************************

  c1->cd(3);
  ntuple->SetMarkerColor(4);
  ntuple->Draw("kineticEnergy");
  gPad->SetLogy();

  //*********************************************************************
  // CREATION OF CSV FILE IN GRAS STYLE
  //*********************************************************************

  Double_t PDGCode;
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t xMom;
  Double_t yMom;
  Double_t zMom;
  Double_t kineticEnergy;

  if ( ! rowWise ) {
    ntuple->SetBranchAddress("PDGCode",&PDGCode);
    ntuple->SetBranchAddress("x",&x);
    ntuple->SetBranchAddress("y",&y);
    ntuple->SetBranchAddress("z",&z);
    ntuple->SetBranchAddress("xMom",&xMom);
    ntuple->SetBranchAddress("yMom",&yMom);
    ntuple->SetBranchAddress("zMom",&zMom);
    ntuple->SetBranchAddress("kineticEnergy",&kineticEnergy);
  }
  else {
    SetLeafAddress(ntuple, "PDGCode",&PDGCode);
    SetLeafAddress(ntuple, "x",&x);
    SetLeafAddress(ntuple, "y",&y);
    SetLeafAddress(ntuple, "z",&z);
    SetLeafAddress(ntuple, "xMom",&xMom);
    SetLeafAddress(ntuple, "yMom",&yMom);
    SetLeafAddress(ntuple, "zMom",&zMom);
    SetLeafAddress(ntuple, "kineticEnergy",&kineticEnergy);
  }

  system ("rm -rf GRAS_tmp.csv");
  std::ofstream csv("GRAS_tmp.csv");
  
  // REMINDER: GRAS CSV FORMAT
  /*
  'eventID','',    1,'Event ID'
  'trackID','',    1,'Track ID'
  'PDGEncoding','',    1,'PDG Encoding of particle type'
  'Z','',    1,'Atomic number of nucleus'
  'A','',    1,'Nucleon number of nucleus'
  'Q','',    1,'Charge'
  'excitation','MeV',    1,'Nuclear excitation energy'
  'weight','',    1,'Particle weight'
  'KE','MeV',    1,'Kinetic energy'
  'T','s',    1,'Time coordinate'
  'XG','mm',    3,'Position coordinate (global)'
  'XL','mm',    3,'Position coordinate (local)'
  'VG','',    3,'Direction cosine (global)'
  'VL','',    3,'Direction cosine (local)'
  'D','mm',    1,'Distance from last volume/scatter'
  'instanceID','',    1,'Instance ID of PV'
  */

  cout << "-> Number of events = " << ntuple->GetEntries() << endl;

  for (Int_t j=0;j<ntuple->GetEntries(); j++)
  {
    ntuple->GetEntry(j);
    
    csv 
    << "       "
    << j+1 << "," << "             "
    << 1 << "," << "             "
    << PDGCode << "," << "             "
    << 1 << "," << "             "
    << 1 << "," << "             "
    << -999 << "," << "             "
    << -999 << "," << "             "
    << -999 << "," << "             "
    << kineticEnergy*1e-3 << ","  << "             "// keV to MeV
    << 0 << "," << "             "
    << x*1e-3 << ","  << "             "// um to mm
    << y*1e-3 << ","  << "             "// um to mm
    << z*1e-3 << ","  << "             "// um to mm
    << -1 << "," << "             "
    << -1 << "," << "             "
    << -1 << "," << "             "
    << xMom << "," << "             "
    << yMom << "," << "             "
    << zMom << "," << "             "
    << -1 << "," << "             "
    << -1 << "," << "             "
    << -1 << "," << "             "
    << 0 << "," << "             "
    << 1 << endl;
  }

  system ("rm -rf GRAS.csv");
  system ("cat GRAS_header.txt GRAS_tmp.csv GRAS_eof.txt >> GRAS.csv");
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
