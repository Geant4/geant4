{  
#include <vector>
    gROOT->Reset();
    ifstream in;
    TFile * file = new TFile("Dose.root","RECREATE");
    // LOAD THE EXPERIMENTAL DATA FILE
    // CONTAINED IN THE DIRECTORY
    // hadrontherapy/experimentalData/proton/BraggPeak


    TNtuple *ntupleExperimental = new TNtuple("ntupleExperimental","Protons, exp. data", "depthExp:EdepExp");

    vector <Float_t> vec_dose, vec_iX;

    TString doseFileExp = "../../../experimentalData/proton/BraggPeak/62MeVInWater.out"; 
    
    cout << "Reading file \" " << doseFileExp << "\" ... ";
    Long64_t nlines = ntupleExperimental -> ReadFile(doseFileExp, "depthExp:EdepExp"); 
    if (nlines <=0){cout << "Error: Check file \"" << doseFileExp << "\"\n"; return;}

    printf("%d Experimental points found\n", nlines); 

    Float_t depthExp, EdepExp;
    ntupleExperimental -> SetBranchAddress("EdepExp", &EdepExp);
    ntupleExperimental -> SetBranchAddress("depthExp", &depthExp);


    // CREATION AND NORMALISATION TO THE FIRST POINT  OF AN NTUPLE CONTAINING THE EXPERIMENTAL DATA
    Int_t nentries = (Int_t)ntupleExperimental -> GetEntries();   
    ntupleExperimental -> GetEntry(0);
    Float_t normFactor = EdepExp;
    for (Int_t l = 0; l<nentries; l++)
    {
	ntupleExperimental -> GetEntry(l);
	vec_dose.push_back(EdepExp);
	vec_iX.push_back(depthExp);
    }

    ntupleExperimental->Reset(); 

    for (Int_t l=0;l<vec_dose.size();l++)
    {
	depthExp = vec_iX[l];
	EdepExp = vec_dose[l]/normFactor;
	ntupleExperimental -> Fill(depthExp, EdepExp);
    }

    //*****************************************************************************
    // Load Simulation file  
    TString doseFileSim = "../../../SimulationOutputs/proton/BraggPeak/Dose.out"; 
    TNtuple *TNtupleSim = new TNtuple("SimTree","dose from ascii file", "iX:jY:kZ:dose"); 

    in.open(doseFileSim);
    if (!in.is_open()){cout << "Error: Check file \"" << doseFileSim << "\"\n"; return;}
    Char_t n[5];
    Float_t f1, f2, f3, f4;
    nlines = 0;
    cout << "Reading file \" " << doseFileSim << "\" ... ";
    // Skip j,j,k,Dose strings
    in >> n >> n >> n >> n;
    do{
	in >> f1 >> f2 >> f3 >> f4;
	nlines++;
	TNtupleSim -> Fill(f1, f2, f3, f4);
	nlines++;}
    while(in.good());

    if (nlines <= 0){cout << "\nNo data found! Check file \"" << doseFileSim << "\"\n"; return;}
    in.close();

    Float_t iX, dose, sumDose = 0., norm = 0. ;
    TNtupleSim -> SetBranchAddress("dose", &dose);
    TNtupleSim -> SetBranchAddress("iX", &iX);

    // Normalize data to 1 at the entry!

    nentries = (Int_t)TNtupleSim -> GetEntries();   
    TNtupleSim -> GetEntry(0);
    Int_t oldX = iX;
    vec_iX.clear();
    vec_dose.clear();
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
    printf("%d Simulated points found\n", vec_iX.size()); 
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
