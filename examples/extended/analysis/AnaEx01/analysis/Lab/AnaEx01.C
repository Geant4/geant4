//
//  To do some analysis over the AnaEx01.root file with the ROOT tool.
//
//  Note that .root is not the default file format and that someone has 
// to ask for an export of a ROOT file in the AnaEx01AnalysisManager.cc 
// constructor.
//
//  Usage :
//    UNIX> root
//    root[0] .X AnaEx01.C
//

{
  TFile* file = new TFile("AnaEx01.root");
  file->ls();

  TTree* tree = (TTree*)file->Get("/tuples/AnaEx01");
  tree->Print();

  TCanvas* canvas = new TCanvas("canvas","EAbs",10,10,800,600);

  tree->Draw("EAbs","","");                                 

  canvas->Update();
}
