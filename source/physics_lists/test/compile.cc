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
#include "G4RunManager.hh" 

//9.4beta #include "FTFC.hh"
//9.4beta #include "FTFP.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_TRV.hh"
//#include "FTFP_BERT_EMV.hh"
//#include "FTFP_BERT_EMX.hh"

#include "LBE.hh"
//9.4beta #include "LHEP_BERT_HP.hh"
//9.4beta #include "LHEP_BERT.hh"
//#include "LHEP_EMV.hh"
//#include "LHEP.hh"
//#include "LHEP_PRECO_HP.hh"

#include "QBBC.hh"

//#include "QGSC_EMV.hh"
//#include "QGSC.hh"
//#include "QGSC_EFLOW.hh"
#include "QGSC_BERT.hh"
//#include "QGSC_QGSC.hh"
//#include "QGSC_CHIPS.hh"
//#include "CHIPS.hh"

//#include "QGSP.hh"
//#include "QGSP_DIF.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BERT.hh"
//#include "QGSP_BERT_EMX.hh"
//#include "QGSP_BERT_EMV.hh"
//#include "QGSP_BERT_CHIPS.hh"
//#include "QGSP_BERT_NOLEP.hh"
//#include "QGSP_BERT_TRV.hh"
//#include "QGSP_BERT_DIF.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC.hh"
//#include "QGSP_BIC_EMY.hh"
#include "G4IonQMDPhysics.hh"
//#include "QGSP_EMV.hh"
//#include "QGSP_EMX.hh"
//#include "QGSP_QEL.hh"
//#include "QGSP_CASC.hh"
//#include "QGSP_INCL_ABLA.hh"
#include "QGSP_INCLXX.hh"

#include "QGSP_FTFP_BERT.hh"

#include "QGS_BIC.hh"
#include "FTF_BIC.hh"

#include "Shielding.hh"

#include "CLHEP/Random/RanluxEngine.h" 

int main(int argc,char** argv) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4RunManager* runManager = new G4RunManager; 
//  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); 

//9.4beta   G4VModularPhysicsList * theFTF1 = new FTFC; 
//9.4beta   G4VModularPhysicsList * theFTF2 = new FTFP;
  G4VModularPhysicsList * theFTF3 = new FTFP_BERT;
  G4VModularPhysicsList * theFTF3a = new FTFP_BERT_TRV;
  //G4VModularPhysicsList * theFTF3b = new FTFP_BERT_EMV;
  //G4VModularPhysicsList * theFTF3c = new FTFP_BERT_EMX;

  G4VModularPhysicsList *theLBE = new LBE; 

//9.4beta   G4VModularPhysicsList *thePL4 = new LHEP_BERT_HP; 
//9.4beta   G4VModularPhysicsList *thePL5 = new LHEP_BERT; 
  //G4VModularPhysicsList *thePL8 = new LHEP_EMV; 
  //G4VModularPhysicsList *thePL10 = new LHEP; 

  G4VModularPhysicsList *thePL15 = new QBBC;

//9.4beta   G4VModularPhysicsList *thePL16 = new QGSC_EMV;
//9.4beta   G4VModularPhysicsList *thePL17 = new QGSC;
//replaced  G4VModularPhysicsList *thePL17a = new QGSC_EFLOW;
  //G4VModularPhysicsList *thePL18 = new QGSC_BERT;
//9.4beta   G4VModularPhysicsList *thePL18a = new QGSC_QGSC;
  //G4VModularPhysicsList *thePL18b = new QGSC_CHIPS;
  //G4VModularPhysicsList *thePL19 = new CHIPS;
    
  G4VModularPhysicsList *thePL20 = new QGSP_BERT_HP; 
  G4VModularPhysicsList *thePL21 = new QGSP_BERT;
  //G4VModularPhysicsList *thePL21a = new QGSP_BERT_EMX;
  //G4VModularPhysicsList *thePL21b = new QGSP_BERT_EMV;
  //G4VModularPhysicsList *thePL21c = new QGSP_BERT_CHIPS;
  //G4VModularPhysicsList *thePL21d = new QGSP_BERT_NOLEP;
  //G4VModularPhysicsList *thePL21e = new QGSP_BERT_TRV;
//9.4beta   G4VModularPhysicsList *thePL21b = new QGSP_BERT_DIF;
  G4VModularPhysicsList *thePL22 = new QGSP_BIC_HP; 
  G4VModularPhysicsList *thePL23 = new QGSP_BIC;
  G4VModularPhysicsList *thePL23a= new QGSP_BIC(1);
      G4IonQMDPhysics * ionQMD=new G4IonQMDPhysics("IonQMD",1);
      thePL23a->RegisterPhysics(ionQMD);
      ionQMD->ConstructProcess();
  //G4VModularPhysicsList *thePL30b = new QGSP_BIC_EMY;
//9.4beta   G4VModularPhysicsList *thePL24 = new QGSP_EMV;
//replaced  G4VModularPhysicsList *thePL25 = new QGSP_EMX;
  //G4VModularPhysicsList *thePL27 = new QGSP;
//9.4beta   G4VModularPhysicsList *thePL27a = new QGSP_DIF;
  //G4VModularPhysicsList *thePL28 = new QGSP;
//  G4VModularPhysicsList *thePL29 = new QGSP_CASC;
//  G4VModularPhysicsList *thePL29 = new QGSP_INCL_ABLA;
  G4VModularPhysicsList *thePL29 = new QGSP_INCLXX;
  
  G4VModularPhysicsList *thePL30 = new QGS_BIC;
  G4VModularPhysicsList *thePL31 = new FTF_BIC;
  
  G4VModularPhysicsList *thePL40 = new QGSP_FTFP_BERT;

  G4VModularPhysicsList *thePL50 = new Shielding;

  delete runManager; 
  return 0; 
} 
