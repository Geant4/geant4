///////////////////////////////////////////////////////////////////////////////
// File: CCalAnalysis.cc
// Description: CCalAnalysis interfaces all user analysis code
///////////////////////////////////////////////////////////////////////////////

#include "G4RunManager.hh" 

#include "CCalAnalysis.hh"
#include "utils.hh"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotter.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#include <AIDA/IManagedObject.h>

//#define debug

CCalAnalysis* CCalAnalysis::instance = 0;
 
CCalAnalysis::CCalAnalysis() :analysisFactory(0), tree(0), tuple(0), 
  energy(0), profile(0) {

  int i=0;
  for (i=0; i<28; i++) {
    hcalE[i] = 0;
    lateralProfile[i] = 0;
  }
  for (i=0; i<49; i++) {ecalE[i] = 0;}
  for (i=0; i<numberOfTimeSlices; i++) {timeHist[i] = 0;}

  analysisFactory = AIDA_createAnalysisFactory();
  if (analysisFactory) {

    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if (treeFactory) {
      // Tree in memory :
      // Create a "tree" associated to an hbook
      const char* opFileptr = getenv("OSCAR_FILENAME");
      G4String opFilestr = "hcaltb96.his";
      if (opFileptr) opFilestr = opFileptr;
      cout << "********************************************" << endl
	   << "* o/p file on " << opFilestr << endl
	   << "********************************************" << endl << endl;
      tree = treeFactory->create(opFilestr, false, false,"hbook");
      if (tree) {
	// Get a tuple factory :
	ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
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

	  tuple = tupleFactory->create("tuple","Event info", tag);
	  
	  assert(tuple);
	  
	  delete tupleFactory;
	}


	IHistogramFactory* histoFactory	= 
	  analysisFactory->createHistogramFactory(*tree);  
	if (histoFactory) {
	  // Create histos :
	  char id[4], ntupletag[50];
	  //Energy deposit in Hcal layers and profiles
	  for (int i = 0; i<28; i++) {
	    sprintf(id, "%d",i+100);
	    sprintf(ntupletag, "Energy Deposit in Hcal Layer%d",i);
	    hcalE[i] = histoFactory->create1D(id, ntupletag, 300, 0.,
					      300000.);
	    sprintf(id, "%d",i+150);
	    sprintf(ntupletag, "Energy profile in Hcal Layer%d",i);
	    lateralProfile[i] = histoFactory->create1D(id, ntupletag, 
						       300, 0., 300000.);
	  }
	  //Energy deposits in Ecal towers
	  for (i = 0; i<49; i++) {
	    sprintf(id, "%d",i+200);
	    sprintf(ntupletag, "Energy Deposit in Ecal Tower%d",i);
	    ecalE[i] = histoFactory->create1D(id, ntupletag, 300, 0., 
					      300000.);
	  }
	  //Total energy deposit
 	  energy  =  histoFactory->create1D("4000", "Total energy deposited",
					    300, 0., 300000.);
	  //Profile
	  profile = histoFactory->create1D("2000", "Profile", 200, 0., 200.);

	  //Time slices	  
	  for (i=0; i<numberOfTimeSlices; i++){
	    sprintf(id, "%d",i+300);
	    sprintf(ntupletag, "Time slice %d",i);
	    timeHist[i] =  histoFactory->create1D(id, ntupletag, 100, 0.,
						  200000.);
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
  if (tree) 
    delete tree;
  if (analysisFactory) 
    delete analysisFactory; // Will delete tree and histos.
}             

CCalAnalysis* CCalAnalysis::getInstance() {
  if (instance == 0) instance = new CCalAnalysis();
  return instance;
}


// This function fill the 1d histogram of the energies in HCal layers
void CCalAnalysis::InsertEnergyHcal(float* v) {
  for (int i=0; i<28; i++) {
    if (hcalE[i]) {
      double x = v[i];
      hcalE[i]->fill(x);
#ifdef debug
      cout << "Fill Hcal histo " << i << " with " << x << endl;
#endif
    }
  }
}


// This function fill the 1d histogram of the energies in ECal layers
void CCalAnalysis::InsertEnergyEcal(float* v) {
  for (int i=0; i<49; i++) {
    if (ecalE[i]) {
      double x = v[i];
      ecalE[i]->fill(x);
#ifdef debug
      cout << "Fill Ecal histo " << i << " with " << x << endl;
#endif
    }
  }
}


// This function fill the 1d histogram of the lateral profile
void CCalAnalysis::InsertLateralProfile(float* v) {
  for (int i=0; i<28; i++) {
    if (lateralProfile[i]) {
      double x = v[i];
      lateralProfile[i]->fill(x);
#ifdef debug
      cout << "Fill Profile Hcal histo " << i << " with " << x << endl;
#endif
    }
  }
}


// This function fill the 1d histogram of the energy 
void CCalAnalysis::InsertEnergy(float v) {
  if (energy) {
    double x = v;
    energy->fill(x);
#ifdef debug
    cout << "Fill Total energy Hcal histo with " << x << endl;
#endif
  }
}


// This function fill the 1d histograms of time profiles
void CCalAnalysis::InsertTime(float* v) {
  for (int j=0; j<numberOfTimeSlices; j++) {
    if (timeHist[j]) {
      double x = v[j];
      timeHist[j]->fill(x);
#ifdef debug
      cout << "Fill Time profile histo " << j << " with " << x << endl;
#endif
    }
  }
}


void CCalAnalysis::setNtuple(float* hcalE, float* ecalE, float elab, 
			     float x, float y, float z, float edep, 
			     float edec, float edhc) {

  ITuple * ntuple = dynamic_cast<ITuple *> ( tree->find("tuple") );
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

  int i;
  ofstream       oFile;
  oFile.open("time.dat");
  for (i=0; i<numberOfTimeSlices; i++){
    if (timeHist[i]) {
      cout << "time slice " << i
	   << " Mean = "    << timeHist[i]->mean()
	   << " sigma = "   << timeHist[i]->rms()
	   << endl;
      oFile<< i  << "  " 
	   << timeHist[i]->mean()  << "  " 
	   << timeHist[i]->rms() << endl;
    }
  }
  oFile.close(); 

  oFile.open("profile.dat");
 
  float w[28] = {2.,2.,
                 3.,3.,3.,3.,3.,3.,
                 6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,6.,
                 8.,8.,8.,8.,8.,8.};

  double cu = 0, edep=0.;
  for (i=0; i<28; i++) {
    if (hcalE[i]) {
      cout << "Histo  " << i
	   << " Mean = "       << hcalE[i]->mean()
	   << " sigma = "      << hcalE[i]->rms()
	   << endl;
      oFile << i  << "  " 
	    << hcalE[i]->mean()          << "  " 
	    << hcalE[i]->rms()           << "  "
	    << lateralProfile[i]->mean() << "  " 
	    << lateralProfile[i]->rms()  << "  " 
	    << endl;
      
      cu += w[i];
      edep += hcalE[i]->mean();
      double wt = hcalE[i]->mean()/w[i];
      if (profile) profile->fill(cu,wt);
    }
  }
  oFile.close(); 
  
  if (tree) tree->commit();

}


// This member is called at the end of every event 
void CCalAnalysis::EndOfEvent(G4int flag) {

  // The plotter is updated only if there is some
  // hits in the event
  if (!flag) return;
}








