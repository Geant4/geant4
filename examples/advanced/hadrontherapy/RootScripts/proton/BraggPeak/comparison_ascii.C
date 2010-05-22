{  
#include <vector>
    gROOT->Reset();

    // LOAD THE EXPERIMENTAL DATA FILE
    // CONTAINED IN THE DIRECTORY
    // hadrontherapy/experimentalData/proton/BraggPeak
    TFile *experimentalFile = new TFile("../../../experimentalData/proton/BraggPeak/62MeVInWater.root","READ");

    // HERE THE ROOT FILE IS INTERPRETED AS A TREE
    TTree *experimentalTree = (TTree*)experimentalFile -> Get("Experimental62MeVInWater");

    TNtuple *ntupleExperimental = new TNtuple("ntupleExperimental","Protons, exp. data", "depthExp:EdepExp");

    Float_t depthExp, EdepExp;
    experimentalTree -> SetBranchAddress("EdepExp", &EdepExp);
    experimentalTree -> SetBranchAddress("depthExp", &depthExp);


    // CREATION AND NORMALISATION TO THE FIRST POINT  OF AN NTUPLE CONTAINING THE EXPERIMENTAL DATA
    Int_t nentries = (Int_t)experimentalTree -> GetEntries();   
    experimentalTree -> GetEntry(0);
    Float_t normFactor = EdepExp;
    for (Int_t l = 0; l<nentries; l++)
    {
	experimentalTree -> GetEntry(l);
	// Is there a method to directly modify data in TTree?
	ntupleExperimental->Fill(depthExp,EdepExp/normFactor); 
    }

//*****************************************************************************
    // Load Simulation file  
    TString doseFile = "../../../SimulationOutputs/proton/BraggPeak/Dose.out"; 
    TFile * file = new TFile("Dose.root","RECREATE");
    TNtuple *TNtupleSim = new TNtuple("SimTree","dose from ascii file", "iX:jY:kZ:dose"); 

       
   ifstream in(doseFile);
   if (!in.is_open()){cout << "Check file \"" << doseFile << "\"\n"; return;}
    Float_t f1,f2,f3,f4;
    Int_t nlines = 0;
    Char_t n[5];
    // Skip j,j,k,Dose strings
    in >> n >> n >> n >> n;
    do{
	in >> f1 >> f2 >> f3 >> f4;
	nlines++;
	TNtupleSim -> Fill(f1, f2, f3, f4);
	nlines++;}
	while(in.good());


    if (nlines <= 0){cout << "No data found! Check file \"" << doseFile << "\"\n"; return;}
    printf("%d points found\n", nlines); 

    Float_t iX, dose, sumDose = 0., norm = 0. ;
    TNtupleSim -> SetBranchAddress("dose", &dose);
    TNtupleSim -> SetBranchAddress("iX", &iX);

    // Normalize data to 1 at the entry!
    
    Int_t nentries = (Int_t)TNtupleSim -> GetEntries();   
    TNtupleSim -> GetEntry(0);
    std::vector <Float_t> vec_dose, vec_iX;
    Int_t oldX = iX;
    // Sum dose along X --> i
    for (Int_t l = 0; l<nentries; l++)
    {
	TNtupleSim -> GetEntry(l);
	if (iX==oldX){ sumDose+=dose;}
	else
	{
	    vec_dose.push_back(sumDose);
	    vec_iX.push_back(oldX);
	    sumDose = dose;
	    oldX = iX;
	}
    }
    // Mean over the first points
    for (Int_t l=0; l<5; l++)
    {
	norm +=vec_dose[l];
    }
    norm /=l;

    TNtupleSim -> Reset();
    // Fill with normalized values. I suppose that slabs/voxel are 0.2 mm depth
    for (Int_t l=0;l<vec_dose.size();l++)
    {
	// Slabs (voxels) are 0.2 mm depth
	iX = 0.1 + (vec_iX[l]*0.2);
	dose = vec_dose[l]/norm;
	TNtupleSim -> Fill(iX, 0, 0, dose);
    }

    TCanvas *c1 = new TCanvas ("c1","c1",200,10,600,400);
    // Triangle 
    TNtupleSim-> SetMarkerStyle(26);
    TNtupleSim-> SetMarkerSize(0.8);
    // Square
    //TNtupleSim-> SetMarkerStyle(25);
    // Star
    //TNtupleSim-> SetMarkerStyle(3);
    // circle
    ntupleExperimental -> SetMarkerStyle(4);
    ntupleExperimental -> SetMarkerColor(2);
    ntupleExperimental -> SetMarkerSize(0.8);

    ntupleExperimental  -> Draw("EdepExp:depthExp");
    TNtupleSim ->  Draw("dose:iX","","same");


    // LEGEND
    leg = new TLegend(0.50,0.60,0.20,0.70); 
    leg -> SetTextSize(0.035);
    leg -> SetFillColor(0);
    leg -> AddEntry(ntupleExperimental, "Experiment", "P");
    leg -> AddEntry(TNtupleSim, "Simulation", "P");
    leg -> Draw();

    
    //c1->SaveAs("braggPeakComparison.pdf");
}
