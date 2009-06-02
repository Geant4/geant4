#include "G4RunManager.hh" 

#include "FTFC.hh"
#include "FTFP.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_TRV.hh"

#include "LBE.hh"
#include "LHEP_BERT_HP.hh"
#include "LHEP_BERT.hh"
#include "LHEP_EMV.hh"
#include "LHEP.hh"
#include "LHEP_PRECO_HP.hh"

#include "QBBC.hh"

#include "QGSC_EMV.hh"
#include "QGSC.hh"
#include "QGSC_BERT.hh"

#include "QGSP.hh"
#include "QGSP_DIF.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_DIF.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_EMV.hh"
#include "QGSP_EMX.hh"
#include "QGSP_QEL.hh"
//#include "QGSP_CASC.hh"

#include "QGS_BIC.hh"
#include "FTF_BIC.hh"

#include "CLHEP/Random/RanluxEngine.h" 

int main(int argc,char** argv) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4RunManager* runManager = new G4RunManager; 
//  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); 

  G4VModularPhysicsList * theFTF1 = new FTFC; 
  G4VModularPhysicsList * theFTF2 = new FTFP;
  G4VModularPhysicsList * theFTF3 = new FTFP_BERT;
  G4VModularPhysicsList * theFTF3a = new FTFP_BERT_TRV;

  G4VModularPhysicsList *theLBE = new LBE; 

  G4VModularPhysicsList *thePL4 = new LHEP_BERT_HP; 
  G4VModularPhysicsList *thePL5 = new LHEP_BERT; 
  G4VModularPhysicsList *thePL8 = new LHEP_EMV; 
  G4VModularPhysicsList *thePL10 = new LHEP; 

  G4VModularPhysicsList *thePL15 = new QBBC;

  G4VModularPhysicsList *thePL16 = new QGSC_EMV;
  G4VModularPhysicsList *thePL17 = new QGSC;
  G4VModularPhysicsList *thePL18 = new QGSC_BERT;

  G4VModularPhysicsList *thePL20 = new QGSP_BERT_HP; 
  G4VModularPhysicsList *thePL21 = new QGSP_BERT;
  G4VModularPhysicsList *thePL21a = new QGSP_BERT_DIF;
  G4VModularPhysicsList *thePL22 = new QGSP_BIC_HP; 
  G4VModularPhysicsList *thePL23 = new QGSP_BIC;
  G4VModularPhysicsList *thePL24 = new QGSP_EMV;
  G4VModularPhysicsList *thePL25 = new QGSP_EMX;
  G4VModularPhysicsList *thePL27 = new QGSP;
  G4VModularPhysicsList *thePL27a = new QGSP_DIF;
  G4VModularPhysicsList *thePL28 = new QGSP_QEL;
//  G4VModularPhysicsList *thePL29 = new QGSP_CASC;
  
  G4VModularPhysicsList *thePL30 = new QGS_BIC;
  G4VModularPhysicsList *thePL31 = new FTF_BIC;
  


  delete runManager; 
  return 0; 
} 
