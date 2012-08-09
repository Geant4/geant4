#include <iostream>
using namespace std;
#include "TMath.h"

void xmacro(TString f1="thin.root", TString f2="thick.root", TString outf="result.pdf")
{
  //Reduce verbosity: Info level is disabled or, since it will send to stderr will
  //create problems...
  gErrorIgnoreLevel=kError;
  gStyle->SetFrameBorderMode(0);  //no frame border
  gStyle->SetCanvasBorderMode(0);  //no canvas border
  gStyle->SetOptTitle(0);  //no title
  	
	
//***********************************************************************************************************************
//  Setup data for Geant4 version 
//***********************************************************************************************************************	

	Double_t sigma_xL94(0.), pTL94(0.);
	Double_t sigma_xS94(0.), pTS94(0.);

	// Thick variables
	TFile *fileL94 = new TFile(f2);
	TTree *TL94 = dynamic_cast<TTree*>(fileL94->Get("sigma;"));
	
	TBranch *SXL94 = TL94->GetBranch("sigma_x");
	//TLeaf *lL93 = bL94->GetLeaf("sigma_z");
	TL94->SetBranchAddress("sigma_x",&sigma_xL94);


	TBranch *b2L94 = TL94->GetBranch("pT");
	//TLeaf *l2L94 = b2L94->GetLeaf("pT");
	TL94->SetBranchAddress("pT",&pTL94);


	// Thin variables
	TFile *fileS94 = new TFile(f1);
	TTree *TS94 = dynamic_cast<TTree*>(fileS94->Get("sigma;"));
	
	TBranch *SXS94 = TS94->GetBranch("sigma_x");
	//TLeaf *lS94 = bS94->GetLeaf("sigma_z");
	TS94->SetBranchAddress("sigma_x",&sigma_xS94);
	TBranch *b2S94 = TS94->GetBranch("pT");
	//TLeaf *l2S94 = b2->GetLeaf("pT");
	TS94->SetBranchAddress("pT",&pTS94);


//***********************************************************************************************************************
//  Setup data for Geant4 version 9.2.p04
//***********************************************************************************************************************
/*
	Double_t sigma_xL92(0.),  pTL92(0.);
	Double_t sigma_xS92(0.),  pTS92(0.);

	// Thick variables
	TFile *fileL92 = new TFile("Thick9.2.root");
	TTree *TL92 = dynamic_cast<TTree*>(fileL92->Get("sigma;"));
	
	TBranch *SXL92 = TL92->GetBranch("sigma_x");
	//TLeaf *lL92 = bL92->GetLeaf("sigma_z");
	TL92->SetBranchAddress("sigma_x",&sigma_xL92);

	TBranch *b2L92 = TL92->GetBranch("pT");
	//TLeaf *l2L92 = b2L92->GetLeaf("pT");
	TL92->SetBranchAddress("pT",&pTL92);


	// Thin variables
	TFile *fileS92 = new TFile("Thin9.2.root");
	TTree *TS92 = dynamic_cast<TTree*>(fileS92->Get("sigma;"));
	
	TBranch *SXS92 = TS92->GetBranch("sigma_x");
	//TLeaf *lS94 = bS94->GetLeaf("sigma_z");
	TS92->SetBranchAddress("sigma_x",&sigma_xS92);

	TBranch *b2S92 = TS92->GetBranch("pT");
	//TLeaf *l2S92 = b2->GetLeaf("pT");
	TS92->SetBranchAddress("pT",&pTS92);

*/
//***********************************************************************************************************************
// Check that writing of events was OK and match
//***********************************************************************************************************************
	

	// Check that the number of entries in each root file match
	Int_t nevent = TL94->GetEntries();
	//Int_t nevent2 = TL92->GetEntries();

	if (nevent != TS94->GetEntries())
	{
		printf("Number of entries in  do not match check the inputloop.sh for each ");
	}
	/*else if (nevent2 != TS92->GetEntries())
	{
		printf("Number of entries in 9.2 do not match check the inputloop.sh for each ");;
	}*/

	/*if(nevent == nevent2)
	{
		printf("No: of Entries:  " << nevent);;
	}
	else
	{
		printf("Different number of Events fo each geant4 version make sure inputloop.sh have same number of input pT: ");;
	}*/


//***********************************************************************************************************************
// Create plots from the data
//***********************************************************************************************************************


	//gStyle->SetOptFit();
	//  Initialise a pointer to the canvas
	TCanvas *c1 = new TCanvas("c1","Various Geant4 IP_x vs 1/pT",200,10,1090,600);
	//TCanvas *c2 = new TCanvas("c2","Various Geant4 IP_y vs 1/pT",200,10,700,500);
	//TCanvas *c3 = new TCanvas("c3","Various Geant4 IP_z vs 1/pT",200,10,700,500);
	//TCanvas *c4 = new TCanvas("c1","Magnitude of IP vs 1/pT",200,10,700,500);
	c1->SetGrid();	
//	c1->Divide(4);	

	TMultiGraph *mgx = new TMultiGraph("Various G4 IPs", "Various G4 IPs");
	//TMultiGraph *mgmag = new TMultiGraph("Various G4 IPs plot of <IP>_{3D} vs 1/p_T", "Various G4 IPs");
        
	/*TH1F *dtxL94 = new TH1F("x", "L94 h1 title", 100, 0, 4.4);
	TH1F *dtxS94 = new TH1F("x", "S94 h1 title", 100, 0, 4.4);
	TH1F *dtxL92 = new TH1F("x", "L92 h1 title", 100, 0, 4.4);
	TH1F *dtxS92 = new TH1F("x", "S92 h1 title", 100, 0, 4.4);*/
	
	TGraph *dtxL94 = new TGraph(nevent);
        TGraph *dtxS94 = new TGraph(nevent);
	/*TGraph *dtxL92 = new TGraph(nevent2);
        TGraph *dtxS92 = new TGraph(nevent2);*/
	
	// Configure style for the graphs version 
	dtxL94->SetMarkerColor(kRed);
	dtxL94->SetMarkerStyle(20);
	dtxL94->SetTitle("Thick ");	
	dtxS94->SetMarkerColor(kBlack);
	dtxS94->SetMarkerStyle(20);
	dtxS94->SetTitle("Thin ");
	

	// Configure style for the graphs version 9.2
	/*dtxL92->SetMarkerColor(kGreen);
	dtxL92->SetMarkerStyle(20);
	dtxL92->SetTitle("Thick v9.2.p04");
	dtxS92->SetMarkerColor(kBlue);
	dtxS92->SetMarkerStyle(20);
	dtxS92->SetTitle("Thin v9.2.p04");*/

	// Define the fitting functions
	TF1 *fitLx94 = new TF1("fitresultxL94", "pol1");
	TF1 *fitSx94 = new TF1("fitresultxS94", "pol1");
	/*TF1 *fitLx92 = new TF1("fitresultxL92", "pol1");
	TF1 *fitSx92 = new TF1("fitresultxS92", "pol1");*/
	

	fitLx94->SetRange(0.0,3.5);
	fitLx94->SetLineColor(kRed);
	fitLx94->SetLineWidth(1);
	fitSx94->SetRange(0.0,3.5);
	fitSx94->SetLineColor(kBlack);
	fitSx94->SetLineWidth(1);
	/*fitLx92->SetRange(0.0,3.5);
	fitLx92->SetLineColor(kGreen);
	fitLx92->SetLineWidth(1);
	fitSx92->SetRange(0.0,3.5);
	fitSx92->SetLineColor(kBlue);
	fitSx92->SetLineWidth(1);*/


		
	for (Int_t i(0); i<nevent; ++i)
	{

		// simga_x and pT entries
		SXL94->GetEvent(i);
		b2L94->GetEvent(i);
		SXS94->GetEvent(i);
		b2S94->GetEvent(i);
		/*SXL92->GetEvent(i);
		b2L92->GetEvent(i);
		SXS92->GetEvent(i);
		b2S92->GetEvent(i);*/
		// sigma_x plots in microns
		dtxL94->SetPoint(i,1/pTL94,1000*sigma_xL94);
		dtxS94->SetPoint(i,1/pTS94,1000*sigma_xS94);
		//dtxL92->SetPoint(i,1/pTL92,1000*sigma_xL92);
		//dtxS92->SetPoint(i,1/pTS92,1000*sigma_xS92);

	}


	// Apply the fits

	/*	dtLmag94 = Fit(fitmagL94, RQM);
		dtSmag94 = Fit(fitmagL94, RQM);
		dtLmag92 = Fit(fitmagL92, RQM);
		dtSmag92 = Fit(fitmagL92, RQM);*/

	// Add TGraph objects to TMultigraph fo multiple plot	

	c1->Clear();
	//c2->Clear();
	//c3->Clear();
	//c4->Clear();
	// Draw graphs

	
	//dtxL94 = Fit(fitLx94, "RQM");
	/*dtxS94 = Fit(fitLx94, "RQM");
	dtxL92 = Fit(fitLx92, "RQM");
	dtxS92 = Fit(fitLx92, "RQM");*/
	//c1->cd(1);
	mgx->Add(dtxL94);
	mgx->Add(dtxS94);
	//mgx->Add(dtxL92);
	//mgx->Add(dtxS92);
	mgx->Draw("ap");
	
	


	//c1->cd(2);
	//dtxS94->Fit("pol1", "RQ", 0, 4.0);
	//dtxS94->Draw("SAME");

	//c1->cd(3);
	//dtxL92->Fit("pol1", "RQ", 0, 4.0);
	//dtxL92->Draw("SAME");

	//c1->cd(4);
	//dtxS92->Fit("pol1", "RQ", 0, 4.0);
	//dtxS92->Draw("SAME");



	// Add a legend to the graph
	leg = new TLegend(0.6,0.7,0.93,0.93);  //coordinates are fractions
	//of pad dimensions
	leg->AddEntry(dtxL94,"Thick ","p");  // "l" means line
	// use "f" for a box
	// use "p" for a point
	leg->AddEntry(dtxS94,"Thin ","p");
	//leg->AddEntry(dtxL92,"Thick v9.2.p04","p");
	//leg->AddEntry(dtxS92,"Thin v9.2.p04","p");
	leg->SetHeader("Geant4 Version");
	leg->Draw();
	
	dtxL94->Fit(fitLx94, "RQM");
	dtxS94->Fit(fitSx94, "RQM");
	//dtxL92->Fit(fitLx92, "RQM");
	//dtxS92->Fit(fitSx92, "RQM");
	
	Double_t fLp0,fLp0_e,fLp1,fLp1_e;
	Double_t fSp0,fSp0_e,fSp1,fSp1_e;
	fLp0=fitLx94->GetParameter(0);
	fLp0_e=fitLx94->GetParError(0);

	fLp1=fitLx94->GetParameter(1);
	fLp1_e=fitLx94->GetParError(1);
	fSp0=fitSx94->GetParameter(0);
	fSp0_e=fitSx94->GetParError(0);
	fSp1=fitSx94->GetParameter(1);
	fSp1_e=fitSx94->GetParError(1);
	
	//Compute Student's t-test on fit results for both parameters
	Int_t ndf = fitLx94->GetNDF() + fitSx94->GetNDF();
	Double_t t0= ( fLp0 - fSp0 ) / TMath::Sqrt( fLp0_e*fLp0_e + fSp0_e*fLp0_e );
	Double_t pval0 = TMath::StudentI( t0 , ndf );//Integral from -inf,t0
	if ( t0 > 0 ) pval0 = 1-pval0;
	Double_t t1= ( fLp1 - fSp1 ) / TMath::Sqrt( fLp1_e*fLp1_e + fSp1_e*fLp1_e );
	Double_t pval1 = TMath::StudentI( t1 , ndf );
	if ( t1 > 0 ) pval1 = 1-pval1;

	cout<<"Student's t-test on intercept: "<<t0<<" p-value="<<pval0<<" (one-tailed)";//<<endl;
	if ( pval0<0.025 ) 
	  {
	    cout<<" test is failing, limit is p-val=0.025";
	    cerr<<"Student's t-test failing for \"Intercept\" fitted value. Check file: "<<outf<<endl;
	  }
	cout<<endl;
	cout<<"Student's t-test on slope: "<<t1<<" p-value="<<pval1<<" (one-tailed)";//<<endl;
	if ( pval1<0.025 ) { 
		    cout<<" test is failing, limit is p-val=0.025";
		    cerr<<"Student's t-test failing for \"Slope\" fitted value. Check file: "<<outf<<endl;
	}
	cout<<endl;

	char textL94[100],textS94[100],textL92[100],textS92[100];
   	sprintf(textL94,"Thick: %g + %g /p_{T}", fitLx94->GetParameter(0), fitLx94->GetParameter(1));
   	sprintf(textS94,"Thin: %g + %g /p_{T}", fitSx94->GetParameter(0), fitSx94->GetParameter(1));
   	//sprintf(textL92," %g + %g /p_{T}", fitLx92->GetParameter(0), fitLx92->GetParameter(1));
   	//sprintf(textS92," %g + %g /p_{T}", fitSx92->GetParameter(0), fitSx92->GetParameter(1));
	
	TPaveText *pt = new TPaveText(0.3,20.5,2.3,25.006);
	//cout<<textL94<<endl;
	//cout<<textS94<<endl;
	//std::cout << textL92 << std::endl;
	//std::cout << textS92 << std::endl;

	pt->AddText(textL94);
	pt->AddText(textS94);
	//pt->AddText(textL92);
	//pt->AddText(textS92);
	pt->Draw();
	
	
	TString Title="IP_x vs 1/p_{T}";
	TString xaxis="1/p_{T} (c/GeV)";
	TString yaxis="#sigma_{IP_{x}}  (10^{4} events)  (#mu m)";
		mgx->SetTitle(Title);
	mgx->GetXaxis()->SetTitle(xaxis);
	mgx->GetXaxis()->SetRangeUser(0.,6.);
	//mgx->GetYaxis()->SetRangeUser(0.,0.0067);
	mgx->GetYaxis()->SetTitle(yaxis);
	mgx->GetXaxis()->SetLabelSize(.03);
	mgx->GetYaxis()->SetLabelSize(.03);
	c1->Update();

	c1->Draw();
	//c1->Draw();
	//c1->Print("G4_IP_Resolution_Study.pdf")

	//***********************************************************************************************************************
	//***********************************************************************************************************************


	//c1->Modified();
        c1->Print(outf);
}

