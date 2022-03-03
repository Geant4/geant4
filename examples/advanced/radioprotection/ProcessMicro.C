{
////////////////////////////////////////////
	gROOT->Reset(); 
	
	bool saveImage = true;
	string imageFileName = "microdosimetricSpectra.png";
	
	bool saveLinearEdep = true;
	bool saveLogEdep = true;
	bool saveydySpec = true;
	
////////////////////////////////////////////	
//Open File where data has been stored
////////////////////////////////////////////
	TFile* f = new TFile("radioprotection_t0.root");

	TDirectory* dir = (TDirectory*)f->Get("radioprotection_ntuple"); 
	TTree * nEdep = (TTree*)dir->Get("102");
	Double_t edep;
	nEdep->SetBranchAddress("edep", &edep);
	
	double meanChordLength = 1; //change based on the thickness of the sensitive volume (make sure to apply a conversion factor if converting to a tissue equivalent material) 
////////////////////////////////////////////
//Put hits in microdosimeter into linear binned space histogram (typical of experimental setup)
////////////////////////////////////////////
	int numberLinearBins = 4096;
	double minLinealEnergy = 0;
	double maxLinealEnergy = 1000;
	TH1D*  histoLin = new TH1D("h1", "f(y)", numberLinearBins, minLinealEnergy, maxLinealEnergy); 
	
	TH1D*  histoYFY = new TH1D("h2", "yf(y)", numberLinearBins, minLinealEnergy, maxLinealEnergy); 
	TH1D*  histoDY = new TH1D("h3", "d(y)", numberLinearBins, minLinealEnergy, maxLinealEnergy);
	
////////////////////////////////////////////
//Create log bin for microdosimetric plotting
////////////////////////////////////////////
	//Log bin stuff 
	int B = 20; //the number of increments in each decade
	//eg if B = 10: 1 2 3 ... 9 10 20 30 ... 90 100 200 ...
	double yMin = 0.01; // y = lineal energy
	double yMax = 1000;

	double range = yMax / yMin;
	double fooRange = range;
	double fooMag = 10; //any old silly no > 1
	int order = 0;
	while (fooMag > 1)
	{
		fooMag = fooRange / 10 ;
		fooRange = fooRange / 10 ;
		order = order + 1 ; 
	}

	int binsPerDecade = B;
	const int numberOfLogBins = binsPerDecade * order +1; //get the total number of bins 

	vector<double> xbins;
	for (int b = 0; b < numberOfLogBins; b++)
	{
		xbins.push_back((yMin)*pow(10, ((b) / double(B))));
		//	cout << b << " " << xbins[b] << endl;
	}

	TH1D*  histoLog = new TH1D("l1", "edepLin", (numberOfLogBins-1), &xbins[0]);


////////////////////////////////////////////
//loop through the hits in the detector stored in the ntuple
////////////////////////////////////////////
	int numberHits = nEdep->GetEntries();
	cout << "Number of hits in detector = " << numberHits << endl;
	
	for (int i = 0; i < numberHits; i++)
	{
		nEdep->GetEntry(i);
		if (edep <= 0)
		{
			cout << "Edep = 0" << endl;
		}
		
		histoLin->Fill(edep/meanChordLength);
		histoLog->Fill(edep/meanChordLength);
	}
	
////////////////////////////////////////////
//Normalise the log histogram based on the bin width
////////////////////////////////////////////	
	for (int bb = 1; bb <= numberOfLogBins-1; bb++)
	{
		histoLog->SetBinContent(bb, (histoLog->GetBinContent(bb) / histoLog->GetBinWidth(bb)));
	}
	
	
////////////////////////////////////////////
//---------Save values to file--------------
////////////////////////////////////////////
	ofstream outFile;
	if (saveLinearEdep)
	{
		outFile.open("linEdepSpec.txt");
		for (int i = 1 ; i <= numberLinearBins; i++)
		{
			outFile << histoLin->GetBinCenter(i) * meanChordLength << "\t" << histoLin->GetBinContent(i) << endl;
		}
		outFile.close();
	}
	if (saveLogEdep)
	{
		outFile.open("logEdepSpec.txt");
		for (int bb = 1; bb <= numberOfLogBins-1; bb++)
		{
			outFile << histoLog->GetBinCenter(bb) * meanChordLength << "\t" << histoLog->GetBinContent(bb) << endl;
		}
		outFile.close();
	}

////////////////////////////////////////////
//----------Normalise Histo to get f(y)----
////////////////////////////////////////////
	double intOfHist = histoLin->Integral("width") ; //gets counts in each histogram and multiplies by bin width
	histoLin->Scale(1./intOfHist);
	//histoLin->Sumw2();
	
	double intOfLogHist = histoLog->Integral("width") ;
	histoLog->Scale(1./intOfLogHist);
	//histoLog->Sumw2();


//////////////////////////////////////////////////////////////
//-------Calculate RBE using the modified MK model-----------
/////////////////////////////////////////////////////////////
//Technically the use of the modified MKM is applied in 
//therapeutic ion beams and not used for radioprotection applications
//though it is included here for convience.
//NOTE: the form used here is for "heavy ions" and not applicable 
//to the more dose depedent protons

	//---------------------------------
	//Constants
	//---------------------------------
	const double satPar = 150.; //saturation paramater obv. in keV/um
	double satParSqr = satPar * satPar;
	const double pi = 3.14159265359;
	const double SiToMuscleCorrection = 0.57;

	//Cell reltated numbers
	double alpha0 = 0.13;//0.164;//0.13; //Gy^-1 //0.13 for carbon 
	double beta = 0.05; //Gy^-2
	double rho = 1.; //g/cm^3
	double rSV = 0.42; //um, the radius of the HSG cell
	//-----------------------------------------------------------
	//In order for consistent units, which in their raw form are:
	//(1/(Gy^2))(keV/g)(cm^3/um^3)
	// (1/Gy^2)(keV/um)(cm^3/g)(1/um^2)
	//------------------------------------------------------------
	//So convert (keV/g) into Gy it is scaled by (1.602E-16/0.001)(J/kg)
	//And to get cm^3 into um^3 just scale by 10E12
	//------------------------------------------------------------
	double scaleGray = 1.602E-16/0.001;
	double scaleLength = 1.E12;

//-------Calculate y*--------
	double topIntegral = 0.;
	double botIntegral = 0.;

	//top integral
	for (Int_t i = 1; i <= numberLinearBins; i++)
	{
		topIntegral += (1 - exp(-(histoLin->GetBinCenter(i) * histoLin->GetBinCenter(i)) / satParSqr))*histoLin->GetBinContent(i) * histoLin->GetBinWidth(i);
		//				(1 - exp(-y^2/y_{0}^2)) * f(y) * dy 
	}
	//bottom integral
	for (Int_t i = 1; i <= numberLinearBins; i++)
	{
		botIntegral += (histoLin->GetBinCenter(i) * histoLin->GetBinContent(i) * histoLin->GetBinWidth(i));
		//				y (size of each being time bin number) * f(y) * dy (bin size)
	}

	//cout << "Top: " << topIntegral << endl;
	//cout << "Bot: " << botIntegral << endl;

	//saturation corrected dose averaged linel energy value (y*)
	double satCor = satParSqr * topIntegral / botIntegral;
	//(y_{0})^2*top/bottom
	//cout << "--------------------------------" << endl;
	//cout << "y* = " << satCor << endl;
	//cout << "--------------------------------" << endl;


//-----------alpha--------
	double alpha = alpha0 + (beta / (rho*pi*rSV*rSV))*satCor*scaleGray*scaleLength;
	//cout << "alpha = " << alpha << endl;

//-----------finally RBE itself--------
	double ln01 = -2.302585093;
	double RBE10;//= (2*beta*5.)((sqrt(alpha*alpha - 4*beta*ln01) - alpha));

	double topRBE = (2 * beta*5.);
	double botRBE = (sqrt(alpha*alpha - 4 * beta*ln01) - alpha);
	RBE10 = topRBE / botRBE;

	//cout << "--------------------------------" << endl;
	cout << "RBE10 = " << RBE10 << endl;
	//cout << "--------------------------------" << endl;
	
////////////////////////////////////////////
//-----------Calculate yF-------------------
////////////////////////////////////////////
//yF = int(y * f(y) dy)
	double linyF = 0.;
	for (int i = 1 ; i <= numberLinearBins; i++)
	{
		linyF += histoLin->GetBinCenter(i) * histoLin->GetBinContent(i) * histoLin->GetBinWidth(i) ;
	}

	cout << "yF = " << linyF << endl;
	
	double logyF = 0;
	for (int i = 1 ; i <= numberOfLogBins; i++)
	{
		logyF += histoLog->GetBinCenter(i) * histoLog->GetBinContent(i) * histoLog->GetBinWidth(i) ;
	}
	cout << "Log yF = " << logyF << endl;


////////////////////////////////////////////
//--------------Calculate d(y)--------------
////////////////////////////////////////////
//d(y) = y*f(y) / yF
	for (int i = 1 ; i <= numberLinearBins; i++)
	{
		histoLin->SetBinContent(i, (histoLin->GetBinCenter(i) * histoLin->GetBinContent(i) / linyF));
	}

	//check that it's normalised to 1
	double intOfHistDy = histoLin->Integral("width") ;
	cout << "Int of d(y) = " << intOfHistDy << endl;
	
	for (int bb = 1; bb <= numberOfLogBins-1; bb++)
	{
		histoLog->SetBinContent(bb, (histoLog->GetBinCenter(bb)* histoLog->GetBinContent(bb) / logyF));
	}
	
	double intOfHistLogDy = histoLog->Integral("width") ;
	cout << "Int of log d(y) = " << intOfHistLogDy << endl;
	
////////////////////////////////////////////
//-------Calculate Average Q(y)-----------
////////////////////////////////////////////
//Q(y) = (a1/y)[1-exp(-(a2)y*y - (a3)y*y*y)]

	//Quality factor values from ICRP 60 report
	double qualityCoeff1 = 5510.;
	double qualityCoeff2 = 5.E-5;
	double qualityCoeff3 = 2.E-7;
	      
	double aveQuality = 0.;
	//AveQuality = int(Q(y)d(y)dy)
	      
	for (int i = 1 ; i <= numberLinearBins; i++)
	{
		aveQuality += (qualityCoeff1/histoLin->GetBinCenter(i))* (1 - exp(-qualityCoeff2*pow(histoLin->GetBinCenter(i), 2.)-qualityCoeff3*pow(histoLin->GetBinCenter(i), 3.)))* histoLin->GetBinContent(i)* histoLin->GetBinWidth(i) ;
	}
	cout << "<Q> = " << aveQuality << endl;

	
////////////////////////////////////////////
//---------Calculate yD--------------------
////////////////////////////////////////////
//yD = int (y * d(y) dy)
	double linyD = 0.;

	for (int i = 1 ; i <= numberLinearBins; i++)
	{
		linyD +=  histoLin->GetBinCenter(i) * histoLin->GetBinContent(i) * histoLin -> GetBinWidth(i) ;
	}
	
	cout << "yD = " << linyD << endl;
	
	double logyD = 0.;
	
	for (int bb = 1; bb <= numberOfLogBins-1; bb++)
	{
		logyD += histoLog->GetBinCenter(bb) * histoLog->GetBinContent(bb) * histoLog -> GetBinWidth(bb) ;
	}
	cout << "log yD = " << logyD << endl;
	
////////////////////////////////////////////
//---------Create yd(y)-------------------
////////////////////////////////////////////
	for (int bb = 1; bb <= numberOfLogBins-1; bb++)
	{
		histoLog->SetBinContent(bb, (histoLog->GetBinCenter(bb)* histoLog->GetBinContent(bb)) );
	}
	//Scale bins by log(10)/B
	//See Appendix B of the ICRU 36 report for more details
	histoLog->Scale((log(10) / B));
	
	if (saveydySpec)
	{
		outFile.open("logYDYspec.txt");
		for (int bb = 1; bb <= numberOfLogBins-1; bb++)
		{
			outFile << histoLog->GetBinCenter(bb) << "\t" << histoLog->GetBinContent(bb) << endl;
		}
		outFile.close();
	}
	
	if (saveImage)
	{
		
		TCanvas *c1 = new TCanvas("c1", "", 1000, 800);
		//gStyle->SetOptStat(0); //un-comment to remove stat box on top right
		
		histoLog->SetTitle("");
		histoLog->GetXaxis()->SetTitle("y (keV/#mum)");
		histoLog->GetYaxis()->SetTitle("yd(y)");
		
		histoLog->GetYaxis()->SetLabelFont(42);
		histoLog->GetXaxis()->SetLabelFont(42);
		histoLog->GetYaxis()->SetTitleFont(42);
		histoLog->GetXaxis()->SetTitleFont(42);
		histoLog->GetYaxis()->CenterTitle(1);
		histoLog->GetXaxis()->CenterTitle(1);
		histoLog->GetXaxis()->SetTitleSize(0.045);
		histoLog->GetYaxis()->SetTitleSize(0.045);
		//histoLog->GetYaxis()->SetRangeUser(0., 1.8);
		//histoLog -> GetXaxis()->SetRangeUser(0.01., 1000.);
		histoLog->GetYaxis()->SetLabelSize(0.04);
		histoLog->GetXaxis()->SetLabelSize(0.04);
		histoLog->SetLineColor(1);
		histoLog->SetLineWidth(4);
		c1->SetLogx();
		//c1->SetLogy();
		//Plotting misbehaves in root6.XX
		histoLog->Draw("");
			
		c1->SaveAs(imageFileName.c_str());
	}
}			
