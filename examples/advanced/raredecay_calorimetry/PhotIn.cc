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
//
// $Id: PhotIn.cc,v 1.10 2006/11/22 09:46:02 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example/advanced/Photon_inefficiency 
//
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

//#define debug

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "PhotInDetectorConstruction.hh"
#include "PhotInDetectorMessenger.hh"
#include "PhotInPhysicsList.hh"
#include "PhotInPrimaryGeneratorAction.hh"
#include "PhotInRunAction.hh"
#include "PhotInEventAction.hh"
#include "PhotInStackingAction.hh"
#include "PhotInSteppingAction.hh"

// Uncomment it if you wish to use old external Reference Physics Lists (+ in GNUmakefile)
//// === Old Reference Physics Lists (external) ===
//#include "FTFC.hh"
//#include "FTFP.hh"
//#include "LBE.hh"
//#include "LHEP.hh"
//#include "LHEP_BERT.hh"
//#include "LHEP_BERT_HP.hh"
//#include "LHEP_BIC.hh"
//#include "LHEP_BIC_HP.hh"
//#include "LHEP_EMV.hh"
////#include "LHEP_GN.hh"
////#include "LHEP_HP.hh"
//#include "LHEP_LEAD.hh"
//#include "LHEP_LEAD_HP.hh"
//#include "LHEP_PRECO.hh"
//#include "LHEP_PRECO_HP.hh"
//#include "QGSC.hh"
//#include "QGSC_LEAD.hh"
//#include "QGSC_LEAD_HP.hh"
//#include "QGSP.hh"
//#include "QGSP.hh"
//#include "QGSP_BERT.hh"
//#include "QGSP_BERT_HP.hh"
//#include "QGSP_BIC.hh"
//#include "QGSP_EMV.hh"
////#include "QGSP_GN.hh"
//#include "QGSP_HP.hh"

// Uncomment it if you wish to use new internal Reference Physics Lists (+ in GNUmakefile)
//#include "G4ModularPhysicsList.hh"
//// === Corresponding New Reference Physics Lists (internal) ===
//#include "G4PL_CHIPS_QGS.hh"
//#include "G4PL_CHIPS_FTF.hh"
//#include "G4PL_CHIPS_LHE.hh"
//#include "G4PL_FTFC_CASP.hh"
//#include "G4PL_FTFC_CASP_GN_HP.hh"
//#include "G4PL_FTFP_CASP.hh"
//#include "G4PL_FTFP_CASP_GN_HP.hh"
//#include "G4PL_UNDERGROUND.hh"
//#include "G4PL_LHEP_CASP.hh"
//#include "G4PL_LHEP_CASP_GN.hh"
//#include "G4PL_LHEP_CASP_HP.hh"
//#include "G4PL_LHEP_CASP_GN_HP.hh"
//#include "G4PL_LHEP_BERT.hh"
//#include "G4PL_LHEP_BERT_HP.hh"
//#include "G4PL_LHEP_BERT_GN_HP.hh"
//#include "G4PL_LHEP_BINC.hh"
//#include "G4PL_LHEP_BINC_HP.hh"
//#include "G4PL_LHEP_BINC_GN_HP.hh"
//#include "G4PL_LHEP_MARS.hh"
//#include "G4PL_LHEP_MARS_HP.hh"
//#include "G4PL_LHEP_MARS_GN_HP.hh"
//#include "G4PL_LHEP_PREC.hh"
//#include "G4PL_LHEP_PREC_HP.hh"
//#include "G4PL_LHEP_PREC_GN_HP.hh"
//#include "G4PL_QGSC_CASP.hh"
//#include "G4PL_QGSC_CASP_GN_HP.hh"
//#include "G4PL_QGSC_MARS.hh"
//#include "G4PL_QGSC_MARS_HP.hh"
//#include "G4PL_QGSC_MARS_GN_HP.hh"
//#include "G4PL_QGSP_CASP.hh"
//#include "G4PL_QGSP_CASP_GN_HP.hh"
//#include "G4PL_QGSP_BERT.hh"
//#include "G4PL_QGSP_BERT_GN_HP.hh"
//#include "G4PL_QGSP_BINC.hh"
//#include "G4PL_QGSP_BINC_GN_HP.hh"
//#include "G4PL_QGSP_CASP_GN.hh"
//#include "G4PL_QGSP_CASP_GN_HP.hh"

int main(int argc,char** argv)
{
#ifdef debug
  G4cout<<"PhotIn: Example is started"<<G4endl;
#endif

 // Construct the default run manager
 G4RunManager* runManager = new G4RunManager; // Del? M.K.

#ifdef debug
  G4cout<<"PhotIn: Run manager is created"<<G4endl;
#endif

  PhotInPhysicsList* physList = new PhotInPhysicsList; // Del? M.K.

#ifdef debug
  G4cout<<"PhotIn: PhotIn Physics List is constructed"<<G4endl;
#endif

  runManager->SetUserInitialization(physList);

#ifdef debug
  G4cout<<"PhotIn: PhotIn Physics List is transfered to Run Manager"<<G4endl;
#endif
  //***LOOKHERE***CHOOSE THE REFERENCE PHYSICS LISTS INSTEAD OF THE PhotInPhysicsList above
  // ==== Initialization of old RPL ====
  //runManager->SetUserInitialization(new LHEP);         // external (ext) LHEP     
  // ==== Initialization of new RPL ====
  //runManager->SetUserInitialization(new G4PL_CHIPS_QGS); // GHIPS+QGSC
  //runManager->SetUserInitialization(new G4PL_CHIPS_FTF); // CHIPS+FTFC
  //runManager->SetUserInitialization(new G4PL_CHIPS_LHE); // CHIPS+LHEP
  //
  //runManager->SetUserInitialization(new G4PL_FTFC_CASP_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new G4PL_FTFC_CASP);// int FTFC_CASP == ext FTFC
  //runManager->SetUserInitialization(new FTFC); // ext FTFC
  //
  //runManager->SetUserInitialization(new G4PL_FTFP_CASP_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new G4PL_FTFP_CASP);// int FTFP_CASP == ext FTFP
  //runManager->SetUserInitialization(new FTFP); // ext FTFP
  //
  //runManager->SetUserInitialization(new G4PL_LHEP_BERT_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new LHEP_BERT_HP); // ext LHEP_BERT_HP
  //runManager->SetUserInitialization(new G4PL_LHEP_BERT);// int LHEP_BERT == ext LHEP_BERT
  //runManager->SetUserInitialization(new LHEP_BERT); // ext LHEP_BERT
  //
  //runManager->SetUserInitialization(new G4PL_LHEP_BINC_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new LHEP_BIC_HP); // ext LHEP_BIC_HP
  //runManager->SetUserInitialization(new G4PL_LHEP_BINC);// int LHEP_BINC == ext LHEP_BIC
  //runManager->SetUserInitialization(new LHEP_BIC); // ext LHEP_BIC
  //
  //runManager->SetUserInitialization(new G4PL_LHEP_CASP_GN_HP);//int _GN_HP (no ext _HP)
  //runManager->SetUserInitialization(new G4PL_LHEP_CASP);// int LHEP_CASP == ext LHEP
  //runManager->SetUserInitialization(new LHEP_EMV); // ext LHEP_EMV (EM_7.0)
  //runManager->SetUserInitialization(new LHEP); // ext LHEP
  //
  //runManager->SetUserInitialization(new G4PL_LHEP_MARS_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new LHEP_LEAD_HP); // ext LHEP_LEAD_HP
  //runManager->SetUserInitialization(new G4PL_LHEP_MARS);// int LHEP_MARS == ext LHEP_LEAD
  //runManager->SetUserInitialization(new LHEP_LEAD); // ext LHEP_LEAD
  //
  //runManager->SetUserInitialization(new G4PL_LHEP_PREC_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new LHEP_PRECO_HP); // ext LHEP_PRECO
  //runManager->SetUserInitialization(new G4PL_LHEP_PREC);//int LHEP_PREC == ext LHEP_PRECO
  //runManager->SetUserInitialization(new LHEP_PRECO); // ext LHEP_PRECO
  //
  //runManager->SetUserInitialization(new G4PL_QGSC_CASP_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new G4PL_QGSC_CASP);// int QGSC_CASP == ext QGSC
  //runManager->SetUserInitialization(new QGSC); // ext QGSC
  //
  //runManager->SetUserInitialization(new G4PL_QGSC_MARS_GN_HP);// int _GN_HP list **BAD**
  //runManager->SetUserInitialization(new QGSC_LEAD_HP); // ext QGSC_LEAD_HP
  //runManager->SetUserInitialization(new G4PL_QGSC_MARS);// int QGSC_MARS == ext QGSC_LEAD
  //runManager->SetUserInitialization(new QGSC_LEAD); // ext QGSC_LEAD
  //
  //runManager->SetUserInitialization(new G4PL_QGSP_BERT_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new QGSP_BERT_HP); // ext QGSP_BERT_HP
  //runManager->SetUserInitialization(new G4PL_QGSP_BERT);// int QGSP_BERT == ext QGSP_BERT
  //runManager->SetUserInitialization(new QGSP_BERT); // ext QGSP_BERT
  //
  //runManager->SetUserInitialization(new G4PL_QGSP_BINC_GN_HP);// int _GN_HP lists
  //runManager->SetUserInitialization(new G4PL_QGSP_BINC);// int QGSP_BINC == ext QGSP_BIC
  //runManager->SetUserInitialization(new QGSP_BIC); // ext QGSP_BIC
  //
  //runManager->SetUserInitialization(new G4PL_QGSP_CASP_GN_HP);// int _GN_HP list
  //runManager->SetUserInitialization(new QGSP_HP); // ext QGSP_HP
  //runManager->SetUserInitialization(new G4PL_QGSP_CASP);// int QGSP_CASP == ext QGSP
  //runManager->SetUserInitialization(new QGSP_EMV); // ext QGSP_EMV (EM_7.0)
  //runManager->SetUserInitialization(new QGSP); // ext QGSP

  // ==== Initialization of Modular Physics List ====
  //runManager->SetUserInitialization(new G4ModularPhysicsList);     // MPL     
  //***endLOOKHERE***

#ifdef debug
  //G4cout<<"PhotIn: Old/New Reference Physics List LHEP is initialized"<<G4endl;
#endif

  // ================= This is a debugging function =====================
  //runManager->SetUserAction(new PhotInSteppingAction);

  // set mandatory initialization classes
		G4double hx=.5*m;
		G4double hy=.5*m;
		G4double hz= m;
  PhotInDetectorConstruction* detector = new PhotInDetectorConstruction(hx,hy,hz); // Del?
  //detector->SetSectionHalfDimensions(0.5*m,0.5*m,m);// This is how it could be changed(-)

#ifdef debug
  G4cout<<"PhotIn: Detector is created with hx="<<hx<<", hy="<<hy<<", hz="<<hz<<G4endl;
#endif
    
  // set user action classes
  PhotInPrimaryGeneratorAction* primaryGenerator = new PhotInPrimaryGeneratorAction;
  runManager->SetUserAction(primaryGenerator);

#ifdef debug
  G4cout<<"PhotIn: Primary Generator is constructed"<<G4endl;
#endif

		//PhotInDetectorMessenger* detectorMessenger=
  //                        new PhotInDetectorMessenger(detector, primaryGenerator);//Del?
  new PhotInDetectorMessenger(detector, primaryGenerator);

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is constructed"<<G4endl;
#endif

  runManager->SetUserInitialization(detector);

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is constructed"<<G4endl;
#endif
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#ifdef debug
  G4cout<<"PhotIn: Visualization Manager is constructed and initialized"<<G4endl;
#endif

#endif

  //PhotInPrimaryGeneratorAction* gen = (PhotInPrimaryGeneratorAction*)
  //                       G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  //if(gen) gen->SetDetector(detector);
  //else G4cout<<"---Warning--- PhotIn: PhotInPrimaryGeneratorAction* ="<<gen<<G4endl;
  // In the example detectorMessenger is created only for DetectorConstruction to change it

  primaryGenerator->SetDetector(detector);

#ifdef debug
  G4cout<<"PhotIn: Detector is transfered to the primary generator"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInRunAction);

#ifdef debug
  G4cout<<"PhotIn: Run is constructed"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInEventAction);

#ifdef debug
  G4cout<<"PhotIn: Event List is constructed"<<G4endl;
#endif

  runManager->SetUserAction(new PhotInStackingAction);
  
#ifdef debug
  G4cout<<"PhotIn: Stacking Actions are constructed"<<G4endl;
#endif

  // Initialize G4 kernel
  // runManager->Initialize(); // (?) (M.K.)
     
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();

#ifdef debug
  G4cout<<"PhotIn: Start UI initialization"<<G4endl;
#endif  

  if (argc==1)   // Define UI session for interactive mode.
  {
    G4UIsession* session=0;
  
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif    

#ifdef debug
    G4cout<<"PhotIn: UITerminal is constructed"<<G4endl;
#endif
      
    UI->ApplyCommand("/control/execute vis.mac");    

#ifdef debug
    G4cout<<"PhotIn: First UI comand /control/execute vis.mac is sent"<<G4endl;
#endif
    session->SessionStart();

#ifdef debug
    G4cout<<"PhotIn: UI session is started"<<G4endl;
#endif

    delete session;

#ifdef debug
    G4cout<<"PhotIn: UI session is finished"<<G4endl;
#endif

  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];

#ifdef debug
    G4cout<<"PhotIn: Batch mode /control/execute "<<argv[1]<<" is started"<<G4endl;
#endif

    UI->ApplyCommand(command+fileName);
  }

  // job termination
#ifdef G4VIS_USE
  delete visManager;

#ifdef debug
  G4cout<<"PhotIn: Visualization Manager is deleted"<<G4endl;
#endif

#endif

  //delete detectorMessenger; // ?

#ifdef debug
  G4cout<<"PhotIn: Detector Messenger is deleted"<<G4endl;
#endif

  delete runManager;

#ifdef debug
  G4cout<<"PhotIn: Run Manager is deleted"<<G4endl;
#endif

  return 0;
}

