{
TFile f("remsim.root");
TDirectory* dir = (TDirectory*)f.Get("histo");
TTree* h2 =(TTree*)dir->Get("1");// Insert the ID of the histo to plot
h2.Draw();
}
