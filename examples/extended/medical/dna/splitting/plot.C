void plot() {
    Int_t numberOfSplit = 200;

    TFile* ref = new TFile("vrt.root");
    TH1F* clusterSize = (TH1F*)ref->Get("2");
    TH1F* energyDeposit = (TH1F*)ref->Get("1");

    // Re-normalization by the number of split to get absolute frequency histogram
    // Not necessary for energy deposit, as we scored edep * weight
    clusterSize->Scale(1.0/numberOfSplit); 

    TCanvas* c1 = new TCanvas("c1","",1200,550);
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetLogy();
    clusterSize->Draw();
    c1->cd(2);
    energyDeposit->Draw();
    c1->Update();
}
