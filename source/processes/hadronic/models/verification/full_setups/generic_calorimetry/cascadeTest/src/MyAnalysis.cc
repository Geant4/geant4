#include "MyAnalysis.hh"
#include <string>
#include "G4RunManager.hh" 

#include <AIDA/AIDA.h>


MyAnalysis* MyAnalysis::instance = 0;

 
MyAnalysis::MyAnalysis() : 
  analysisFactory(0), tree(0), tuple(0) {

  //ALB analysisFactory = AIDA_createAnalysisFactory(); //***LOOKHERE***
  if (analysisFactory) {    
    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if (treeFactory) {
      // Tree in memory: create a "tree" associated to an hbook file.
      bool readOnly = false; // we want to write.
      bool createNew = true; // create file if it doesn't exist.

      tree = treeFactory->create("my.his", "hbook", readOnly, createNew);

      if (tree) {
	// Get a tuple factory :
	AIDA::ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
	if (tupleFactory) {
	  // Create a tuple :
	  std::string tag = "float ID, E, EDEP_ACT, EDEP_CAL";
	  //tuple = tupleFactory->create("1","Event info", tag); // Column wise (default)
	  tuple = tupleFactory->create("1","Event info", tag, "--preferRWN"); // Row wise
	  assert(tuple);
	  delete tupleFactory;
	}
      }
      delete treeFactory; // It will not delete the ITree.
    }
  }
}


MyAnalysis::~MyAnalysis() {
  if (tree) {
    delete tree;
    tree = 0;
  }
  if (analysisFactory) {
    delete analysisFactory;
    analysisFactory = 0;
  }
}


MyAnalysis* MyAnalysis::getInstance() {
  if (instance == 0) instance = new MyAnalysis();
  return instance;
}


void MyAnalysis::init() {
  if (tuple) tuple->reset();
}                       


void MyAnalysis::finish() {             
  if (tree) {
    tree->commit();
    tree->close();
  }
}


void MyAnalysis::fillNtuple( float incidentParticleId, 
			     float incidentParticleEnergy, 
			     float totalEnergyDepositedInActiveLayers,
			     float totalEnergyDepositedInCalorimeter ) {
  if (tuple) {
    tuple->fill( tuple->findColumn( "ID" ), incidentParticleId );
    tuple->fill( tuple->findColumn( "E" ), incidentParticleEnergy );
    tuple->fill( tuple->findColumn( "EDEP_ACT" ), totalEnergyDepositedInActiveLayers );
    tuple->fill( tuple->findColumn( "EDEP_CAL" ), totalEnergyDepositedInCalorimeter );
    tuple->addRow();
  }
}

