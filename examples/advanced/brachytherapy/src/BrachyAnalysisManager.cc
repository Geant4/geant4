

#include <stdlib.h>
#include "g4std/fstream"
#include "BrachyAnalysisManager.hh"

#include "G4ios.hh"

#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"

#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/ITuple.h"

BrachyAnalysisManager* BrachyAnalysisManager::instance = 0;

BrachyAnalysisManager::BrachyAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0)
  

{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  ITreeFactory     * treeFact = aFact->createTreeFactory();
 
 
  
 //parameters for the TreeFactory
  bool fileExists = false;
  bool readOnly   = false;
  std::string fileName="Brachy3.hbk";
  theTree = treeFact->create(fileName, readOnly, fileExists, "hbook");

  delete treeFact;
  //HistoFactory and TupleFactory depend on theTree
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree );
 
}



BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
  delete tupFact;
  tupFact=0;

   delete histFact;
  histFact=0;

  delete theTree;
  histFact=0;

  delete aFact;
  aFact = 0;

}

BrachyAnalysisManager* BrachyAnalysisManager::getInstance()
{
  if (instance == 0) instance = new BrachyAnalysisManager;
  return instance;
}


void BrachyAnalysisManager::book() 
{

  // histograms and ntuple are managed by theTree
  // h2 e h3 are not necessary  ,useful for:check for non-zero
  //you could do  histFact->create2D("20","Energy, pos",300 ,-150.,150.,300,-15  //0.,150.);
 

 IHistogram2D *h2 = histFact->create2D("20","Energy, pos",300 ,-150.,150.,300,-150.,150.);
 
// check for non-zero

 IHistogram1D *h3 = histFact->create1D("30","Initial Energy", 100,0.,1.);
 // check for non-zero

 //Ntuple management
 std::string columnNames = "float energy, float x, float z";
 std::string options = "";
 ITuple *tup = tupFact->create("1","brachy",columnNames, options);
 // check for non-zero ...
 if (tup) G4cout<<"The Ntuple is non-zero"<<G4endl;
}
 

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void BrachyAnalysisManager::analyse(G4double xx,G4double zz,G4float en)
{

  ITuple * ntuple = dynamic_cast<ITuple *> ( theTree->find("1") );

 
  ntuple->fill(1, en);// fill ( int column, double value )
  ntuple->fill(2, xx);
  ntuple->fill(3, zz);

  // alternatively(if all the columns are of the same type):
  // std::vector<float> row;
  // row.push_back(en);
  // row.push_back(xx);
  // row.push_back(zz);
  // ntuple->fill(row);

  // write Ntuple-row to file
  // Should be called after fill is called for the columns.
  // Unfilled columns will be filled with the default value for that column.
  ntuple->addRow();

}

void BrachyAnalysisManager::hist(G4double x,G4double z, G4float enn)
{ 
  

  IHistogram2D* h2 = dynamic_cast<IHistogram2D *> ( theTree->find("20") );
  h2->fill(x,z,enn);  
}


void BrachyAnalysisManager::Spectrum(G4double Init_En)
{ 
  IHistogram1D* h3 = dynamic_cast<IHistogram1D *> ( theTree->find("30") );
  h3->fill(Init_En);
}


void BrachyAnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}












