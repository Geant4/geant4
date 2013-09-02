{

TFile* f = new TFile("out.root", "UPDATE");
gFile->GetListOfFree()->ls();
f->Close();

TFile* f2 = new TFile("out.root", "UPDATE");
gFile->GetListOfFree()->ls();

}
