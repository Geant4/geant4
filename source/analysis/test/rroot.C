
//  This script shows how to read with CERN-ROOT the out.root
// file containing some histos and a TTree produced by
// the wroot.cpp program.

void rroot() {

  TFile* f = new TFile("out.root");
  if(!f) {
    std::cout << "TFile out.root not found." << std::endl;    
    return;
  }

  TCanvas* plotter = new TCanvas("canvas","",10,10,800,600);
  plotter->Divide(2,2);  

  ////////////////////////////////////////////////////////
  /// histos /////////////////////////////////////////////
  ////////////////////////////////////////////////////////
 {TDirectory* dir = (TDirectory*)f->Get("histo");
  if(!dir) {
    std::cout << "TDirectory histos not found." << std::endl;    
    return;
  }
  TH1D* hrg = (TH1D*)dir->Get("rg");
  if(!hrg) {
    std::cout << "TH1D rg not found." << std::endl;    
    return;
  }
  plotter->cd(1);
  hrg->Draw();

  TProfile* hprof = (TProfile*)dir->Get("prof");
  if(!hprof) {
    std::cout << "TProfile prof not found." << std::endl;    
    return;
  }
  plotter->cd(2);
  hprof->Draw();

  TH2D* hrgbw = (TH2D*)dir->Get("rgbw");
  if(!hrgbw) {
    std::cout << "TH2D rgbw not found." << std::endl;    
    return;
  }
  plotter->cd(3);
  hrgbw->Draw();}

  ////////////////////////////////////////////////////////
  /// TTree produced with wroot::ntuple //////////////////
  ////////////////////////////////////////////////////////
 {TH1D* hrg = new TH1D("hrg","Rand gauss",100,-5,5);

  TTree* tree = (TTree*)f->Get("rg_rbw");
  if(!tree) {
    std::cout << "TTree rg_rbw not found." << std::endl;    
    return;
  }

 {TObjArray* brs = tree->GetListOfBranches();
  for(int i=0;i<brs->GetEntries();i++) {
    TBranch* b = (TBranch*)brs->At(i);
    std::cout << "branch : " << b->GetName() << std::endl;    
  }}

  int index;
  tree->SetBranchAddress("index",&index);
  double rgauss;
  tree->SetBranchAddress("rgauss",&rgauss);
  float rbw;
  tree->SetBranchAddress("rbw",&rbw);

  int num = tree->GetEntries();  
  std::cout << "number of events = " << num << std::endl;

  for(int i=0;i<num;i++) {
    tree->GetEntry(i);
    if(index!=i) {
      std::cout << "problem to read index branch at entry " << i
                << ". Found index was " << index << "."
                << std::endl;
      break;
    }
    //std::cout << "debug : " << rgauss << std::endl;
    hrg->Fill(rgauss,1);
  }

  plotter->cd(4);
  hrg->Draw();}

  plotter->Update();

}
