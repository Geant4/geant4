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
//
// $Id: FCALAnalysisManager.cc
// GEANT4 tag $Name:
//
// Author: Alex Howard (a.s.howard@ic.ac.uk)
//
// History:
// -----------
// 3 Oct 2002 Alex Howard     Created
// 5 Nov 2002 Alex Howard    Successfully Modified to AIDA 3.x
//
// -------------------------------------------------------------------
#ifdef  G4ANALYSIS_USE

#include "FCALAnalysisManager.hh"

//#include "g4std/iomanip"

FCALAnalysisManager* FCALAnalysisManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager::FCALAnalysisManager() :
  af(0), tree(0), hf(0), tpf(0), pf(0)
{
  // tree is created and booked inside book()
  ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager::~FCALAnalysisManager() 
{

  delete pf;
  pf=0;

  delete tpf;
  tpf=0;

  delete hf;
  hf=0;

  delete tree;
  tree=0;

  delete af;
  af = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisManager* FCALAnalysisManager::getInstance()
{
  if (instance == 0) instance = new FCALAnalysisManager;
  return instance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void FCALAnalysisManager::book()

{
  //  histoManager->selectStore("FCAL.his");
  G4String histogramfile = "FCAL.his";

  G4cout << " Histogramfile: " << histogramfile << G4endl;


  //build up  the  factories
  af = AIDA_createAnalysisFactory();


 //parameters for the TreeFactory
  G4bool fileExists = false;
  G4bool readOnly   = false;

  AIDA::ITreeFactory     * tf = af->createTreeFactory();

  tree = tf->create(histogramfile, "hbook", readOnly, fileExists);

  G4cout << "Tree store : " << tree->storeName() << G4endl;

  G4cout << " Booked Hbook File " << G4endl;

  //HistoFactory and TupleFactory depend on theTree
  hf = af->createHistogramFactory( *tree );
  tpf  = af->createTupleFactory(*tree );

 // ---- primary ntuple ------

  AIDA::ITuple* ntuple1 = tpf->create( "1", "Number out of World", 
			     "float OutOfWorld,i,j" );

  assert(ntuple1);

  // ---- secondary ntuple ------   

  AIDA::ITuple* ntuple2 = tpf->create( "2", "Secondary Info", 
				 "float Secondary,i,j");

  assert(ntuple2);

  // ---- tertiary ntuple ------   

 AIDA::ITuple* ntuple3 = tpf->create( "3", "Energy Deposits", 
				"float EmEdep,HadEdep" );

 assert(ntuple3);

  // Creating an 1-dimensional histogram in the root directory of the tree

  AIDA::IHistogram1D* hOutOfWorld;
  hOutOfWorld    = hf->createHistogram1D("10","Number Of OutOfWorld",  100,0.,100.);

  AIDA::IHistogram1D* hSecondary;
  hSecondary      = hf->createHistogram1D("20","Number Of Secondaries", 100,0.,100.);
  
  AIDA::IHistogram1D* hEmEdep;
  hEmEdep  = hf->createHistogram1D("30","Electromagnetic Energy /MeV", 100,0.,100.);
  
  AIDA::IHistogram1D* hHadEdep;
  hHadEdep  = hf->createHistogram1D("30","Hadronic Energy /MeV", 100,0.,100.);
  
  delete tf;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::finish()
{

  // Committing the transaction with the tree
  G4std::cout << "Committing..." << G4std::endl;
  // write all histograms to file
  tree->commit();

  G4std::cout << "Closing the tree..." << G4std::endl;

  // close (will again commit)
  tree->close();

  // extra delete as objects are created in book() method rather than during
  // initialisation of class
  delete pf;
  delete tpf;
  delete hf;
  delete tree;
  delete af;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void FCALAnalysisManager::NumOutOfWorld(G4double OutOfWorld, G4int i, G4int j)

{
  AIDA::IHistogram1D* h1 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("10") );
  h1->fill(OutOfWorld);  // fill(x,y,weight)     
  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("1") );

 // Fill the ntuple
  ntuple->fill( ntuple->findColumn( "OutOfWorld"   ), (G4float) OutOfWorld  );
  ntuple->fill( ntuple->findColumn( "i"  ),           (G4float) i    );
  ntuple->fill( ntuple->findColumn( "j"   ),      (G4float) j   );
//ntuple->fill( ntuple->findColumn( "Event" ),  static_cast<float>(event_id) );

  //Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::Secondaries(G4double Secondary, G4int i, G4int j)

{
  AIDA::IHistogram1D* h2 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("20") );
  h2->fill(Secondary);  // fill(x,y,weight)     
  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("2") );

 // Fill the ntuple
  ntuple->fill( ntuple->findColumn( "Secondary"   ), (G4float) Secondary  );
  ntuple->fill( ntuple->findColumn( "i"  ),           (G4float) i    );
  ntuple->fill( ntuple->findColumn( "j"   ),      (G4float) j   );
//ntuple->fill( ntuple->findColumn( "Event" ),  static_cast<float>(event_id) );

  //Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisManager::Edep(G4double EmEdep, G4double HadEdep)

{
  //EM:
  AIDA::IHistogram1D* h3 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("30") );
  h3->fill(EmEdep);  // fill(x,y,weight)     
  //Had:
  AIDA::IHistogram1D* h4 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("40") );
  h4->fill(HadEdep);  // fill(x,y,weight)     
  AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("3") );

 // Fill the ntuple
  ntuple->fill( ntuple->findColumn( "EmEdep"   ), (G4float) EmEdep  );
  ntuple->fill( ntuple->findColumn( "HadEdep"  ), (G4float) HadEdep );
//ntuple->fill( ntuple->findColumn( "Event" ),  static_cast<float>(event_id) );

  //Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif







