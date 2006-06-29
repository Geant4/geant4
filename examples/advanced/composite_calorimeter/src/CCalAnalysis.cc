//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalAnalysis.cc
// Description: CCalAnalysis interfaces all user analysis code
///////////////////////////////////////////////////////////////////////////////

#ifdef G4ANALYSIS_USE

#include "G4RunManager.hh" 

#include "CCalAnalysis.hh"
#include "CCalutils.hh"

#include <AIDA/AIDA.h>

#include <typeinfo>

//#define debug

CCalAnalysis* CCalAnalysis::instance = 0;
 
CCalAnalysis::CCalAnalysis() :analysisFactory(0), tree(0), tuple(0), energy(0) {

  for (int i=0; i<28; i++) {hcalE[i] = 0;}
  for (int i=0; i<70; i++) {lateralProfile[i] = 0;}
  for (int i=0; i<49; i++) {ecalE[i] = 0;}
  for (int i=0; i<numberOfTimeSlices; i++) {timeHist[i] = 0;}
  for (int i=0; i<2; i++)  {timeProfile[i] = 0;}

  analysisFactory = AIDA_createAnalysisFactory();
  if (analysisFactory) {

    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if (treeFactory) {
      // Tree in memory :
      // Create a "tree" associated to an hbook
      const char* opFileptr = getenv("CCAL_FILENAME");
      G4String opFilestr = "ccal.his";
      if (opFileptr) opFilestr = opFileptr;
      G4cout << "********************************************" << G4endl
	     << "* o/p file on " << opFilestr << G4endl
	     << "********************************************" << G4endl 
	     << G4endl;
      bool readOnly = false; // we want to write.
      bool createNew = true; // create file if it doesn't exist.
      tree = treeFactory->create(opFilestr, "hbook", readOnly,createNew);
      if (tree) {
	// Get a tuple factory :
	AIDA::ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
	if (tupleFactory) {
	  // Create a tuple :
	  G4String tag2, tag = "float";
	  for (int i=0; i<28; i++) {
	    tag2 = tag + " hcal" + i + ",";
	    tag  = tag2;
	  }
	  for (int i=0; i<49; i++) {
	    tag2 = tag + " ecal" + i + ",";
	    tag  = tag2;
	  }
	  tag2 = tag  + " ELAB, XPOS, YPOS, ZPOS";
	  tag  = tag2 + ", EDEP, EDEC, EDHC";

	  tuple = tupleFactory->create("1","Event info", tag); // Column wise (default)
	  //tuple = tupleFactory->create("1","Event info", tag, "--preferRWN"); // Row wise
	  
	  assert(tuple);
	  
	  delete tupleFactory;
	}

	AIDA::IHistogramFactory* histoFactory	= 
	  analysisFactory->createHistogramFactory(*tree);  
	if (histoFactory) {
	  // Create histos :
	  char id[4], ntupletag[50];
	  //Energy deposit in Hcal layers 
	  for (int i = 0; i<28; i++) {
	    sprintf(id, "%d",i+100);
	    sprintf(ntupletag, "Energy Deposit in Hcal Layer%d   in GeV",i);
	    hcalE[i] = histoFactory->createHistogram1D(id, ntupletag, 100, 0., 1.0);
	  }
	  // Energy deposits in Ecal towers
	  for (int i = 0; i<49; i++) {
	    sprintf(id, "%d",i+200);
	    sprintf(ntupletag, "Energy Deposit in Ecal Tower%d   in GeV",i);
	    ecalE[i] = histoFactory->createHistogram1D(id, ntupletag, 100, 0., 100.0);
	  }
	  // Total energy deposit
 	  energy  =  histoFactory->createHistogram1D("4000", "Total energy deposited   in GeV", 
						     100, 0., 100.0);

	  // Time slices	  
	  for (int i=0; i<numberOfTimeSlices; i++){
	    sprintf(id, "%d",i+300);
	    sprintf(ntupletag, "Time slice %d nsec energy profile   in GeV",i);
	    timeHist[i] =  histoFactory->createHistogram1D(id, ntupletag, 100, 0., 100.0);
	  }

	  // Profile of lateral energy deposit in Hcal
	  for (int i = 0; i<70; i++) {
	    sprintf(id, "%d",i+500);
	    sprintf(ntupletag, "Lateral energy profile at %d cm  in GeV",i);
	    lateralProfile[i] = histoFactory->createHistogram1D(id, ntupletag, 100, 0., 10.0);
	  }

	  // Time profile 
	  timeProfile[0] = histoFactory->createHistogram1D("901", "Time Profile in Sensitive Detector", 200, 0., 200.);
	  timeProfile[1] = histoFactory->createHistogram1D("902", "Time Profile in Sensitive+Passive", 200, 0., 200.);

	  delete histoFactory;
	}
    
      }
      delete treeFactory; // Will not delete the ITree.
    }

  }

}


CCalAnalysis::~CCalAnalysis() {
  Finish();
}


void CCalAnalysis::Init() {
}                       

void CCalAnalysis::Finish() {
  if (tree) { 
    delete tree;
    tree = 0;
  }
  if (analysisFactory) {
    delete analysisFactory; // Will delete tree and histos.
    analysisFactory = 0;
  }
}             

CCalAnalysis* CCalAnalysis::getInstance() {
  if (instance == 0) instance = new CCalAnalysis();
  return instance;
}


// This function fill the 1d histogram of the energies in HCal layers
void CCalAnalysis::InsertEnergyHcal(float* v) {
#ifdef debug
  double totalFilledEnergyHcal = 0.0;
#endif      
  for (int i=0; i<28; i++) {
    if (hcalE[i]) {
      double x = v[i];
      hcalE[i]->fill(x);
#ifdef debug
      G4cout << "Fill Hcal histo " << i << " with " << x << G4endl;
      totalFilledEnergyHcal += x;
#endif      
    }    
  }
#ifdef debug
  G4cout << "CCalAnalysis::InsertEnergyHcal: Total filled Energy Hcal histo " 
	 << totalFilledEnergyHcal << G4endl;
#endif      
}


// This function fill the 1d histogram of the energies in ECal layers
void CCalAnalysis::InsertEnergyEcal(float* v) {
#ifdef debug
  double totalFilledEnergyEcal = 0.0;
#endif      
  for (int i=0; i<49; i++) {
    if (ecalE[i]) {
      double x = v[i];
      ecalE[i]->fill(x);
#ifdef debug
      G4cout << "Fill Ecal histo " << i << " with " << x << G4endl;
      totalFilledEnergyEcal += x;
#endif
    }
  }
#ifdef debug
  G4cout << "CCalAnalysis::InsertEnergyEcal: Total filled Energy Ecal histo " 
	 << totalFilledEnergyEcal << G4endl;
#endif      
}


// This function fill the 1d histogram of the lateral profile
void CCalAnalysis::InsertLateralProfile(float* v) {
#ifdef debug
  double totalFilledProfileHcal = 0.0;
#endif
  for (int i=0; i<70; i++) {
    if (lateralProfile[i]) {
      double x = v[i];
      lateralProfile[i]->fill(x);
#ifdef debug
      G4cout << "Fill Profile Hcal histo " << i << " with " << x << G4endl;
      totalFilledProfileHcal += x;
#endif
    }
  }
#ifdef debug
  G4cout << "CCalAnalysis::InsertLateralProfile: Total filled Profile Hcal"
	 << " histo " << totalFilledProfileHcal << G4endl;
#endif      
}


// This function fill the 1d histogram of the energy 
void CCalAnalysis::InsertEnergy(float v) {
  if (energy) {
    double x = v;
    energy->fill(x);
#ifdef debug
    G4cout << "CCalAnalysis::InsertEnergy: Fill Total energy Hcal histo with " 
	   << x << G4endl;
#endif
  }
}


// This function fill the 1d histograms of time profiles at stepping action
void CCalAnalysis::InsertTime(float* v) {
#ifdef debug
  double totalFilledTimeProfile = 0.0;
#endif
  for (int j=0; j<numberOfTimeSlices; j++) {
    if (timeHist[j]) {
      double x = v[j];
      timeHist[j]->fill(x);
#ifdef debug
      G4cout << "Fill Time slice histo " << j << " with " << x << G4endl;
      totalFilledTimeProfile += x;
#endif
    }
    if (timeProfile[1]) {
      double x = v[j];
      double t = j + 0.5;
      timeProfile[1]->fill(t,x);
#ifdef debug
      G4cout << "Fill Time profile histo 1 with " << t << " " << x << G4endl;
#endif
    }
  }
#ifdef debug
  G4cout << "CCalAnalysis::InsertTime: Total filled Time profile histo " 
	 << totalFilledTimeProfile << G4endl;
#endif      
}


// This function fill the 1d histograms of time profiles in SD
void CCalAnalysis::InsertTimeProfile(int hit, double time, double edep) {
  
  if (timeProfile[0]) {
    timeProfile[0]->fill(time,edep);
#ifdef debug
    G4cout << "CCalAnalysis:: Fill Time Profile with Hit " << hit
	   << " Edeposit " << edep << " Gev at " << time << " ns" << G4endl;
#else
    hit=0;  // Just to avoid compiler warning!
#endif
  }
}


void CCalAnalysis::setNtuple(float* hcalE, float* ecalE, float elab, 
			     float x, float y, float z, float edep, 
			     float edec, float edhc) {

  if (tuple) {
    char tag[10];
    for (int i=0; i<28; i++) {
      sprintf (tag, "hcal%d", i);
      tuple->fill(tuple->findColumn(tag),hcalE[i]);
    }
    for (int i=0; i<49; i++) {
      sprintf (tag, "ecal%d", i);
      tuple->fill(tuple->findColumn(tag),ecalE[i]);
    }
    tuple->fill(tuple->findColumn("ELAB"),elab);
    tuple->fill(tuple->findColumn("XPOS"),x);
    tuple->fill(tuple->findColumn("YPOS"),y);
    tuple->fill(tuple->findColumn("ZPOS"),z);
    tuple->fill(tuple->findColumn("EDEP"),edep);
    tuple->fill(tuple->findColumn("EDEC"),edec);
    tuple->fill(tuple->findColumn("EDHC"),edhc);
    tuple->addRow();
#ifdef debug
    G4cout << "CCalAnalysis:: Fill Ntuple " << G4endl;
#endif
  }
}


/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void CCalAnalysis::BeginOfRun(G4int )  { 
  
  /*
  if (energy) energy->reset();
  for (int i=0; i<28; i++) {
    if (hcalE[i]) hcalE[i]->reset();
  }
  for (int i=0; i<70; i++) {
    if (lateralProfile[i]) lateralProfile[i]->reset();
  }
  for (int i=0; i<49; i++) {
    if (ecalE[i]) ecalE[i]->reset();
  }
  */
  for (int i=0; i<numberOfTimeSlices; i++) {
    if (timeHist[i]) timeHist[i]->reset();
  }

}


//  This member is called at the end of each run 
void CCalAnalysis::EndOfRun(G4int )  {

  if (tree) {
    tree->commit();
    tree->close();
  }

}


// This member is called at the end of every event 
void CCalAnalysis::EndOfEvent(G4int flag) {

#ifdef debug
  G4cout << " Check if empty histograms " << G4endl;
  if ( energy ) {
    if ( energy->allEntries() == 0 ) {
      G4cout << "EMPTY HISTO  energy " << G4endl;
    } else if ( energy->allEntries() == energy->extraEntries() ) {
      G4cout << "EXTRA entries only HISTO  energy " << G4endl;
    }
  } else {
      G4cout << "UNDEFINED HISTO  energy " << G4endl;
  }
  for (int i=0; i<28; i++) {
    if ( hcalE[i] ) {
      if ( hcalE[i]->allEntries() == 0 ) {
	G4cout << "EMPTY HISTO  hcal " << i << G4endl;
      } else if ( hcalE[i]->allEntries() == hcalE[i]->extraEntries() ) {
	G4cout << "EXTRA entries only HISTO  hcal " << i << G4endl;
      }
    } else {
      G4cout << "UNDEFINED HISTO  hcal " << i << G4endl;
    }
  }
  for (int i=0; i<70; i++) {
    if ( lateralProfile[i] ) {
      if ( lateralProfile[i]->allEntries() == 0 ) {
	G4cout << "EMPTY HISTO  lateralProfile " << i << G4endl;
      } else if ( lateralProfile[i]->allEntries() == 
		  lateralProfile[i]->extraEntries() ) {
	G4cout << "EXTRA entries only HISTO  lateralProfile " << i << G4endl;
      }
    } else {
      G4cout << "UNDEFINED HISTO  lateralProfile " << i << G4endl;
    }
  }
  for (int i=0; i<49; i++) {
    if ( ecalE[i] ) {
      if ( ecalE[i]->allEntries() == 0 ) {
	G4cout << "EMPTY HISTO  ecalE " << i << G4endl;
      } else if   ( ecalE[i]->allEntries() == ecalE[i]->extraEntries() ) {
	G4cout << "EXTRA entries only HISTO  ecal " << i << G4endl;
      }
    } else {
      G4cout << "UNDEFINED HISTO  hcal " << i << G4endl;
    }
  }
  for (int i=0; i<numberOfTimeSlices; i++) {
    if ( timeHist[i] ) {
      if ( timeHist[i]->allEntries() == 0 ) {
	G4cout << "EMPTY HISTO  timeHist " << i << G4endl;
      } else if ( timeHist[i]->allEntries() == timeHist[i]->extraEntries() ) {
	G4cout << "EXTRA entries only HISTO  timeHist " << i << G4endl;
      }
    } else {
      G4cout << "UNDEFINED HISTO  timeHist " << i << G4endl;
    }
  }
#endif  

  // The plotter is updated only if there is some
  // hits in the event
  if (!flag) return;
}

#endif
