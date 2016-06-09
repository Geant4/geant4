//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

//#define debug

CCalAnalysis* CCalAnalysis::instance = 0;
 
CCalAnalysis::CCalAnalysis() :analysisFactory(0), tree(0), tuple(0), energy(0) {

  int i=0;
  for (i=0; i<28; i++) {
    hcalE[i] = 0;
  }
  for (i=0; i<70; i++) {
    lateralProfile[i] = 0;
  }
  for (i=0; i<49; i++) {ecalE[i] = 0;}
  for (i=0; i<numberOfTimeSlices; i++) {timeHist[i] = 0;}

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
	   << "********************************************" << G4endl << G4endl;
      bool readOnly = false; // we want to write.
      bool createNew = true; // create file if it doesn't exist.
      tree = treeFactory->create(opFilestr, "hbook", readOnly,createNew);
      if (tree) {
	// Get a tuple factory :
	AIDA::ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
	if (tupleFactory) {
	  // Create a tuple :
	  G4String tag2, tag = "float";
	  for (i=0; i<28; i++) {
	    tag2 = tag + " hcal" + i + ",";
	    tag  = tag2;
	  }
	  for (i=0; i<49; i++) {
	    tag2 = tag + " ecal" + i + ",";
	    tag  = tag2;
	  }
	  tag2 = tag  + " ELAB, XPOS, YPOS, ZPOS";
	  tag  = tag2 + ", EDEP, EDEC, EDHC";

	  //tuple = tupleFactory->create("1","Event info", tag); // Column wise (default)
	  tuple = tupleFactory->create("1","Event info", tag, "--preferRWN"); // Row wise
	  
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
	  for (i = 0; i<49; i++) {
	    sprintf(id, "%d",i+200);
	    sprintf(ntupletag, "Energy Deposit in Ecal Tower%d   in GeV",i);
	    ecalE[i] = histoFactory->createHistogram1D(id, ntupletag, 100, 0., 100.0);
	  }
	  // Total energy deposit
 	  energy  =  histoFactory->createHistogram1D("4000", "Total energy deposited   in GeV", 
					    100, 0., 100.0);

	  // Time slices	  
	  for (i=0; i<numberOfTimeSlices; i++){
	    sprintf(id, "%d",i+300);
	    sprintf(ntupletag, "Time slice %d nsec energy profile   in GeV",i);
	    timeHist[i] =  histoFactory->createHistogram1D(id, ntupletag, 100, 0., 100.0);
	  }

	  // Profile of lateral energy deposit in Hcal
	  for (i = 0; i<70; i++) {
	    sprintf(id, "%d",i+400);
	    sprintf(ntupletag, "Lateral energy profile at %d cm  in GeV",i);
	    lateralProfile[i] = histoFactory->createHistogram1D(id, ntupletag, 100, 0., 10.0);
	  }

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
      G4cout << "\t total filled Energy Hcal histo " << totalFilledEnergyHcal << G4endl;
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
      G4cout << "\t total filled Energy Ecal histo " << totalFilledEnergyEcal << G4endl;
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
  G4cout << "\t total filled Profile Hcal histo " << totalFilledProfileHcal << G4endl;
#endif      
}


// This function fill the 1d histogram of the energy 
void CCalAnalysis::InsertEnergy(float v) {
  if (energy) {
    double x = v;
    energy->fill(x);
#ifdef debug
    G4cout << "Fill Total energy Hcal histo with " << x << G4endl;
#endif
  }
}


// This function fill the 1d histograms of time profiles
void CCalAnalysis::InsertTime(float* v) {
#ifdef debug
  double totalFilledTimeProfile = 0.0;
#endif
  for (int j=0; j<numberOfTimeSlices; j++) {
    if (timeHist[j]) {
      double x = v[j];
      timeHist[j]->fill(x);
#ifdef debug
      G4cout << "Fill Time profile histo " << j << " with " << x << G4endl;
      totalFilledTimeProfile += x;
#endif
    }
  }
#ifdef debug
  G4cout << "\t total filled Time profile histo " << totalFilledTimeProfile << G4endl;
#endif      
}


void CCalAnalysis::setNtuple(float* hcalE, float* ecalE, float elab, 
			     float x, float y, float z, float edep, 
			     float edec, float edhc) {

  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("1") );
  if (ntuple) {
    char tag[10];
    int i=0;
    for (i=0; i<28; i++) {
      sprintf (tag, "hcal%d", i);
      ntuple->fill(tuple->findColumn(tag),hcalE[i]);
    }
    for (i=0; i<49; i++) {
      sprintf (tag, "ecal%d", i);
      ntuple->fill(tuple->findColumn(tag),ecalE[i]);
    }
    ntuple->fill(tuple->findColumn("ELAB"),elab);
    ntuple->fill(tuple->findColumn("XPOS"),x);
    ntuple->fill(tuple->findColumn("YPOS"),y);
    ntuple->fill(tuple->findColumn("ZPOS"),z);
    ntuple->fill(tuple->findColumn("EDEP"),edep);
    ntuple->fill(tuple->findColumn("EDEC"),edec);
    ntuple->fill(tuple->findColumn("EDHC"),edhc);
    ntuple->addRow();
  }
}


/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void CCalAnalysis::BeginOfRun(G4int n)  { 
  
  int i=0;  
  /*
  if (energy) energy->reset();
  for (i=0; i<28; i++) {
    if (hcalE[i]) hcalE[i]->reset();
  }
  for (i=0; i<70; i++) {
    if (lateralProfile[i]) lateralProfile[i]->reset();
  }
  for (i=0; i<49; i++) {
    if (ecalE[i]) ecalE[i]->reset();
  }
  */
  for (i=0; i<numberOfTimeSlices; i++) {
    if (timeHist[i]) timeHist[i]->reset();
  }

}


//  This member is called at the end of each run 
void CCalAnalysis::EndOfRun(G4int n)  {

  if (tree) {
    tree->commit();
    tree->close();
  }

}


// This member is called at the end of every event 
void CCalAnalysis::EndOfEvent(G4int flag) {

#ifdef debug
  G4cout << " Check if empty histograms " << G4endl;
  int i=0;  
  if ( energy ) {
    if ( energy->allEntries() == 0 ) {
      G4cout << "EMPTY HISTO  energy " << G4endl;
    } else if ( energy->allEntries() == energy->extraEntries() ) {
      G4cout << "EXTRA entries only HISTO  energy " << G4endl;
    }
  } else {
      G4cout << "UNDEFINED HISTO  energy " << G4endl;
  }
  for (i=0; i<28; i++) {
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
  for (i=0; i<70; i++) {
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
  for (i=0; i<49; i++) {
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
  for (i=0; i<numberOfTimeSlices; i++) {
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
