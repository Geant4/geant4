void macro( TString dirname = "thin" )
{
  //cout<<"Processing directory:"<<dirname<<endl;
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetFillColor(0);  //white fill color
	gStyle->SetFrameBorderMode(0);  //no frame border
	gStyle->SetCanvasBorderMode(0);  //no canvas border
	gStyle->SetOptTitle(0);  //no title
	//gStyle->SetOptStat(0);  //no statistics _or_
	gStyle->SetOptStat(1111111111);
	// -->all statistics (the more 1's y
	// connect file and tree
	//TFile *file = new TFile(dirname+"IPbinpT_0.4.root");
	//TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(dirname+"IPbinpT_0.33.root");
	//if (f==0) {
	//	std::cout << "What This file doesn't exist MAAAAN!" << std::endl;
	//}
	//TTree *tree = (TTree*)gDirectory->Get("sigma"); 
	//   Create the chain as above
	//
	TChain *chain = new TChain("sigma");
	chain->Add(dirname+"/*.root");

	chain->Merge(dirname+".root");
	//cout<<"Merged file: "<<dirname<<".root"<<endl;
	//exit();
	
	//chain->Draw("bx1:pTi");

	/*	for (Int_t i(0); i<9; ++i)
		{
		chain->GetEvent(i);
		}
	 */
	/*chain.SetBranchAddress("event", &event);

	//   Start main loop on all events
	//   In case you want to read only a few branches, use TChain::SetBranchStatus
	//   to activate/deactivate a branch.
	Int_t nevent = chain.GetEntries();
	for (Int_t i=0;i<nevent;i++) {
	chain.GetEvent(i);              //read complete accepted event in memory
	hnseg->Fill(event->GetNseg());  //Fill histogram with number of segments
	}

	//  Draw the histogram
	hnseg->Draw();*/
}
