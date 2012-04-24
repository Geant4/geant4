// 
#ifndef BEAMTESTRUNACTION_HH
#define BEAMTESTRUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "BeamTestRunActionMessenger.hh"

#ifdef G4ANALYSIS_USEROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "Riostream.h"
#include "TObject.h"
#include "TH1F.h"
#include "TKey.h"
#include "TF1.h"
#include "TStyle.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string.h>
#include <sstream>

class G4Run;

class BeamTestRunAction : public G4UserRunAction 
{

	public:

		// Constructor
		BeamTestRunAction(/*Parameters* parameter*/);

		// Destructor
		virtual ~BeamTestRunAction();

		// Methods  
		virtual G4Run* GenerateRun();

		void BeginOfRunAction(const G4Run*);
		void EndOfRunAction(const G4Run*);

		//const G4double getPT() const {return pT;}
		// Fill events to be written to file
		void FillEvents(G4double, G4double, G4double, G4double);
        void SetpT( G4double pt ) { pT = pt; }
        void SetEnergyDeposit(G4double val );
        void AddSuffix( G4String s ) { fNameSuffix << s; }
        void Initialize();
        void MergeFiles( G4String pattern );
        void SetMergedFilename( G4String fn ) { mergedFileName = fn; }
        void Finalize(G4String macroname);
	private:
        std::stringstream fNameSuffix;
        G4String mergedFileName;
        void Write();
        G4double pT;
#ifdef G4ANALYSIS_USEROOT
		TFile *f;
		TString pTString;

		TH1F *ipPerTrackb_x;
		TH1F *ipPerTrackb_y;
		TH1F *ipPerTrackb_z;
		TH1F *ipPerTrackb_mag;
        TH1F *EDepPerTrack;

		//TTree* newTree;
		TTree* sigma;
		Double_t msigma_x, msigma_y, msigma_z, edep;

		TF1* fitresultx;
		TF1* fitresulty;
		TF1* fitresultz;
#endif
        BeamTestRunActionMessenger* messager;
};

#endif

