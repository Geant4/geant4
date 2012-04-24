// ********************************************************************
//
// Matthew Reid 
//
#ifndef BEAMTESTEVENTACTION_HH
#define BEAMTESTEVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
// Call Root headers

// define Hit of Silicon Monitor
#include "BeamTestSiliconMonitorHit.hh"
//#include "BeamTestParameters.hh"
//#include "BeamTestConversion.hh"
#include "BeamTestRunAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include "G4ThreeVector.hh"

class BeamTestRunAction;

class BeamTestEventAction : public G4UserEventAction {

	public:

		// Constructor
		explicit BeamTestEventAction(/*Parameters* parameter,*/ BeamTestRunAction*);

		// Destructor
		~BeamTestEventAction();

		// Metohds
		void BeginOfEventAction(const G4Event*);
		void EndOfEventAction(const G4Event*);
		G4double ScatteredAngle(G4double x_i, G4double z_FinalExitPos);
        void Initialize(G4ThreeVector mom, G4ThreeVector pos);
        void SetNumberOfChambers( G4int val ) { numberOfChambers = val; }
        void AbortEvent() { abortEvent=true; }
	private:
		
		BeamTestRunAction* runAct;
		// Data member
		G4int fHitsCollectionID;
		G4int fHitsCollectionID_monitor;
		G4double p;
		G4double pT;
		G4double angle;
		G4double pz;
		G4ThreeVector I;
		G4ThreeVector E;
		G4ThreeVector PD;
		G4ThreeVector b;
        G4int numberOfChambers;
        G4bool abortEvent;
};

std::string makeFilename( const std::string& basename, double index );

#endif


