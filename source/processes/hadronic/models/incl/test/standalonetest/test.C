{
	TFile *cpp = new TFile("test.root");
	TFile *fort = new TFile("testfort.root");

	TTree *c = (TTree*) cpp->Get("h101");
	TTree *f = (TTree*) fort->Get("h101");

	TCanvas *c1 = new TCanvas();

	c->Draw("Massini");
	f->SetMarkerStyle(2);
	f->Draw("Massini", "", "P, same");
}
