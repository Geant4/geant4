{

//Open File where data has been stored!
TFile f("radioprotection_NEW.root");
TDirectory* dir = f.Get("radioprotection_ntuple"); 

TNtuple * ntuple1 = (TNtuple*)dir->Get("101");   
TNtuple * ntuple2 = (TNtuple*)dir->Get("102"); 
TNtuple * ntuple3 = (TNtuple*)dir->Get("103"); 

int numberOfBinsX = 500;
int numberOfBinsY = 500;
int numberOfBinsZ = 500;

int Xmin = 0;
int Xmax = 1000;
int Ymin = 0;
int Ymax = 1000;
int Zmin = 0;
int Zmax = 1000;

//the type of histogram, how many variable, min and max values plus size of bins
TH1F* edep1Distribution = new TH1F("h0", "Primary Particle Energy Spectrum; Energy (Mev);Frequency",
				     numberOfBinsX, Xmin, Xmax);       //Edep // binning, xmin, xmax, along x direction

//the type of histogram, how many variable, min and max values plus size of bins
TH1F* edep1Distribution = new TH1F("h1", "Energy deposition; Edep (kev);Frequency",
				     numberOfBinsX, Xmin, Xmax);       //Edep // binning, xmin, xmax, along x direction


TH2F* edep2DDistribution = new TH2F("h2", "Energy deposited by Ions; Z; Edep (keV)", 
				     numberOfBinsX, Xmin, Xmax,    // Z     // binning, xmin, xmax, along x direction
				     numberOfBinsY, Ymin, Ymax); 	//Edep   // binning, xmin, xmax, along y direction


TH3F* edep3DDistribution = new TH3F("h3", "3Dedep; Edep (keV) ; A; Z", 
				     numberOfBinsX, Xmin, Xmin,  // edep  // binning, xmin, xmax, along x direction
				     numberOfBinsY, Ymin, Ymax, 	// A     //binning, xmin, xmax, along y direction
				     numberOfBinsZ, Zmin, Zmax);	// Z     //binning, xmin, xmax, along z direction


//Plot Primary Energy of Incident Particle
ntuple1.Draw("Ek>>h0","","");
//Plot Energy Deposition within SV
//ntuple2.Draw("edep>>h1", "", "");
//Plot 2D/3D Histogram of energy with particle type using A and Z
//ntuple3.Draw("Z:edep>>h2", "", "");
//ntuple3.Draw("Z:A:edep>>h3", "", "");
}			
