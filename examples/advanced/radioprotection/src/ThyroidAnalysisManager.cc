

#include <stdlib.h>
#include "g4std/fstream"
#include "ThyroidAnalysisManager.hh"

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

ThyroidAnalysisManager* ThyroidAnalysisManager::instance = 0;

ThyroidAnalysisManager::ThyroidAnalysisManager() : 
  aFact(0), theTree(0), histFact(0), tupFact(0)
  

{
  //build up  the  factories
  aFact = AIDA_createAnalysisFactory();

  ITreeFactory     * treeFact = aFact->createTreeFactory();
 
 
  
 //parameters for the TreeFactory
  bool fileExists = false;
  bool readOnly   = false;
  std::string fileName="Thyroid.hbk";
  theTree = treeFact->create(fileName, readOnly, fileExists, "hbook");

  delete treeFact;
  //HistoFactory and TupleFactory depend on theTree
  histFact = aFact->createHistogramFactory( *theTree );
  tupFact  = aFact->createTupleFactory    ( *theTree );
 
}



ThyroidAnalysisManager::~ThyroidAnalysisManager() 
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

ThyroidAnalysisManager* ThyroidAnalysisManager::getInstance()
{
  if (instance == 0) instance = new ThyroidAnalysisManager;
  return instance;
}


void ThyroidAnalysisManager::book() 
{

  // histograms and ntuple are managed by theTree
  // h2 e h3 are not necessary  ,useful for:check for non-zero
  //you could do  histFact->create2D("20","Energy, pos",300 ,-150.,150.,300,-15  //0.,150.);
 

 IHistogram2D *h2 = histFact->create2D("20","Energy, pos",200 ,0.,10.,600,0.,30.);
 
// check for non-zero

 IHistogram1D *h3 = histFact->create1D("30","Initial Energy", 500,0.00,0.5);
 // check for non-zero

 //Ntuple management
 std::string columnNames = "float energy, float x,float y, float z";
 std::string options = "";
 ITuple *tup = tupFact->create("1","Thyroid",columnNames, options);
 // check for non-zero ...
 if (tup) G4cout<<"The Ntuple is non-zero"<<G4endl;
}
 

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ThyroidAnalysisManager::analyse(G4double xx,G4double yy,G4double zz,G4float en)
{

  ITuple * ntuple = dynamic_cast<ITuple *> ( theTree->find("1") );

 
  ntuple->fill(1, en);// fill ( int column, double value )
  ntuple->fill(2, xx);
 ntuple->fill(3, yy);
  ntuple->fill(4, zz);

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

void ThyroidAnalysisManager::hist(G4double x,G4double y, G4float enn)
{ 
  

  IHistogram2D* h2 = dynamic_cast<IHistogram2D *> ( theTree->find("20") );
  h2->fill(x,y,enn);  
}


void ThyroidAnalysisManager::Spectrum(G4double Init_En)
{ 
  IHistogram1D* h3 = dynamic_cast<IHistogram1D *> ( theTree->find("30") );
  h3->fill(Init_En);
}


void ThyroidAnalysisManager::finish() 
{  
  // write all histograms to file
  theTree->commit();

  // close (will again commit)
  theTree->close();

}












