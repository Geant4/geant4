#include "G4RunManager.hh" 

#include "FTFC.hh"
#include "FTFP.hh"
#include "LBE.hh"
#include "LHEP_BERT_HP.hh"
#include "LHEP_BERT.hh"
#include "LHEP_EMV.hh"
#include "LHEP.hh"
#include "LHEP_PRECO_HP.hh"

#include "QBBC.hh"

#include "QGSC_EMV.hh"
#include "QGSC.hh"

#include "QGSP.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_EMV.hh"
#include "QGSP_EMX.hh"
#include "QGSP_QEL.hh"

#include "CLHEP/Random/RanluxEngine.h" 

int main(int argc,char** argv) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4RunManager* runManager = new G4RunManager; 
//  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); 

  FTFC *thePL1 = new FTFC; 
  FTFP *thePL2 = new FTFP; 

  LBE *thePL3 = new LBE; 

  LHEP_BERT_HP *thePL4 = new LHEP_BERT_HP; 
  LHEP_BERT *thePL5 = new LHEP_BERT; 
  LHEP_EMV *thePL8 = new LHEP_EMV; 
  LHEP *thePL10 = new LHEP; 

  QBBC *thePL15 = new QBBC;

  QGSC_EMV *thePL16 = new QGSC_EMV;
  QGSC *thePL17 = new QGSC;

  QGSP_BERT_HP *thePL20 = new QGSP_BERT_HP; 
  QGSP_BERT *thePL21 = new QGSP_BERT;
  QGSP_BIC_HP *thePL22 = new QGSP_BIC_HP; 
  QGSP_BIC *thePL23 = new QGSP_BIC;
  QGSP_EMV *thePL24 = new QGSP_EMV;
  QGSP_EMX *thePL25 = new QGSP_EMX;
  QGSP *thePL27 = new QGSP;
  QGSP_QEL *thePL28 = new QGSP_QEL;


  delete runManager; 
  return 0; 
} 
