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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISAAnalysisManager class                                        *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#ifdef G4ANALYSIS_USE

#include "LISAAnalysisManager.hh"
#include "globals.hh"

LISAAnalysisManager* LISAAnalysisManager::instance = 0;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

LISAAnalysisManager::LISAAnalysisManager() :
  af(0),
  tf(0), 
  run_tree(0),
  run_tpf(0),
  run_tuple(0)
{;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

LISAAnalysisManager::~LISAAnalysisManager() {

  if(af) {
    delete tf;
    G4cout << " LISAAnalysis -- deleted tree factory" << G4endl;
    delete af;
    G4cout << " LISAAnalysis -- deleted analysis factory" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISAAnalysisManager* LISAAnalysisManager::getInstance() {

  if (!instance) instance = new LISAAnalysisManager();
  return instance;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAAnalysisManager::Dispose() {

  if(instance) {
    delete instance;
    instance = 0;
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAAnalysisManager::Init() {

  G4cout << G4endl << "******* Analysis with AIDA 3.0 *********" << G4endl;

  // create analysis factory
  if( (af = AIDA_createAnalysisFactory()) )
    G4cout << " LISAAnalysis -- created analysis factory" << G4endl;
  
  // create tree factory
  if( (tf = af->createTreeFactory()) )
    G4cout << " LISAAnalysis -- created tree factory" << G4endl;

}




//
// RUN ANALYSIS **************************************************************
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAAnalysisManager::bookRun(G4String hbookfile) {


  // create RunTree
  G4bool fileExists = false;
  G4bool readOnly   = false;
  run_tree = tf->create(hbookfile, "hbook", readOnly, fileExists);
  G4cout << " LISAAnalysis -- tree store: " << run_tree->storeName()<<G4endl;

  // create TupleFactory
  run_tpf = af->createTupleFactory(*run_tree );
  G4cout << " LISAAnalysis -- created NTuple factory" << G4endl;

  // Run Information
  run_tuple = run_tpf->create( "1", "Run Tuple", 
              "float evt,tm,energy,charge,in,out,seed1,seed2");
  assert(run_tuple);
  G4cout << " LISAAnalysis -- created Run NTuple" << G4endl;
  

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAAnalysisManager::FinishRun() {

  // Committing transaction with the tree
  G4cout << "  LISAAnalysis -- committing run_tree..." << G4endl;
  run_tree->commit();
  run_tree->close();

  delete run_tpf;
  delete run_tree;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISAAnalysisManager::analyseRun(G4int evt, G4int tm, G4double energy, 
 G4int charge, G4int in, G4int out, long seed1, long seed2) {

  AIDA::ITuple* ntuple = dynamic_cast<AIDA::ITuple *> ( run_tree->find("1") );

  // Fill the ntuple
  ntuple->fill( ntuple->findColumn( "evt"    ), (float) evt    );
  ntuple->fill( ntuple->findColumn( "tm"     ), (float) tm     );
  ntuple->fill( ntuple->findColumn( "energy" ), (float) energy );
  ntuple->fill( ntuple->findColumn( "charge" ), (float) charge );
  ntuple->fill( ntuple->findColumn( "in"     ), (float) in     );
  ntuple->fill( ntuple->findColumn( "out"    ), (float) out    );
  ntuple->fill( ntuple->findColumn( "seed1"  ), (float) seed1  );
  ntuple->fill( ntuple->findColumn( "seed2"  ), (float) seed2  );

  // Values of attributes are prepared; store them to the nTuple:
  ntuple->addRow();

}


#endif
