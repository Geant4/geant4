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
// * Neither the authors 262of this software system, nor their employing *
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
// $Id: exampleS5.cc,v 1.1 2006-11-28 15:19:42 tnikitin Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01RunAction.hh"
#include "ExN01EventAction.hh"
#include "ExN01SteppingAction.hh"

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UniformMagField.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4MagIntegratorDriver.hh"  
#include <fstream>
///
///MagField
///

///
///EndMagField
///

///
///SetupStepper
///
G4MagIntegratorStepper* SetupStepper(G4int type,G4Mag_UsualEqRhs* pE)
{
  G4MagIntegratorStepper* pStepper;
  switch ( type )
    {
      case 0: pStepper = new G4ExplicitEuler( pE ); break;
      case 1: pStepper = new G4ImplicitEuler( pE ); break;
      case 2: pStepper = new G4SimpleRunge( pE ); break;
      case 3: pStepper = new G4SimpleHeum( pE ); break;
      case 4: pStepper = new G4ClassicalRK4( pE ); break;
      case 5: pStepper = new G4HelixExplicitEuler( pE ); break;
      case 6: pStepper = new G4HelixImplicitEuler( pE ); break;
      case 7: pStepper = new G4HelixSimpleRunge( pE ); break;
      case 8: pStepper = new G4CashKarpRKF45( pE );    break;
      case 9: pStepper = new G4ExactHelixStepper( pE );    break;
      case 10: pStepper = new G4RKG3_Stepper( pE );    break;
      default: pStepper = 0;
    }
  return pStepper;
}
///
///EndSetupStepper
///

///
///CompareStepStepper
///
void CompareStep(G4double y[],G4double dydx[],G4double hhh,
                 G4double yout1[],G4double yout2[],G4double ydif[],
                 G4double yerr1[],G4double yerr2[],
                  G4MagIntegratorStepper* pSt1,G4MagIntegratorStepper* pSt2)
{   
  //G4double yerr1[G4FieldTrack::ncompSVEC];
  //   G4double yerr2[G4FieldTrack::ncompSVEC];
     G4double ydiferr[G4FieldTrack::ncompSVEC];
     G4double r[10];
     //
     //r0=r;r1=Pr;r2=difR;r3=difPr;r4=diferrR;r5=diferrPr;r6=err1r;r7=err1Pr;
     //r8=err2r;r9=err2Pr;
     //
     r[0]=std::sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
     r[1]=std::sqrt(y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
    pSt1->Stepper(y,dydx,hhh,yout1,yerr1);
    //G4cout << " One step with  stepper type 1"<<G4endl;
    //G4cout << " Step Length= "<<hhh/mm<<"mm"<<G4endl;
    //       G4cout << "Y    DyDx      Yout "<<G4endl;
    //G4cout.precision(10);
    //for (G4int i=0;i<6;i++)
    //G4cout << y[i]<<"  "<<dydx[i]<<"  "<<yout1[i]<<G4endl;
    pSt2->Stepper(y,dydx,hhh,yout2,yerr2);
    //G4cout << " One step with  stepper type 2 "<<G4endl;
    //G4cout << "Y    DyDx      Yout "<<G4endl;
    //G4cout.precision(10);
    //for (G4int i=0;i<6;i++)  
    // G4cout << y[i]<<"  "<<dydx[i]<<"  "<<yout2[i]<<G4endl;
     for (G4int i=0;i<6;i++) ydif[i]= std::abs(yout1[i]-yout2[i]);
     for (G4int i=0;i<6;i++) ydiferr[i]=std::sqrt(yerr1[i]*yerr1[i]+yerr2[i]*yerr2[i]);
       r[2]=std::sqrt(ydif[0]*ydif[0]+ydif[1]*ydif[1]+ydif[2]*ydif[2]);
       r[3]=std::sqrt(ydif[3]*ydif[3]+ydif[4]*ydif[4]+ydif[5]*ydif[5]);
       r[4]=std::sqrt(ydiferr[0]*ydiferr[0]+ydiferr[1]*ydiferr[1]+ydiferr[2]*ydiferr[2]);
       r[5]=std::sqrt(ydiferr[3]*ydiferr[3]+ydiferr[4]*ydiferr[4]+ydiferr[5]*ydiferr[5]);
       r[6]=std::sqrt(yerr1[0]*yerr1[0]+yerr1[1]*yerr1[1]+yerr1[2]*yerr1[2]);
       r[7]=std::sqrt(yerr1[3]*yerr1[3]+yerr1[4]*yerr1[4]+yerr1[5]*yerr1[5]);
       r[8]=std::sqrt(yerr2[0]*yerr2[0]+yerr2[1]*yerr2[1]+yerr2[2]*yerr2[2]);
       r[9]=std::sqrt(yerr2[3]*yerr2[3]+yerr2[4]*yerr2[4]+yerr2[5]*yerr2[5]);
    G4cout.precision(10);
    for (G4int i=0;i<6;i++)  
      // G4cout << "deltaY["<<i<<"]="<<ydif[i]/mm<<"mm"<<" "
      //	     <<ydiferr[i]/mm<<" "
      //        <<yerr1[i]/mm<<" "<<yerr2[i]/mm<<G4endl;
      G4cout <<hhh<<" "<<y[i]<<"  "<<ydif[i]/mm <<" "
     	     <<ydiferr[i]/mm<<" "
            <<yerr1[i]/mm<<" "<<yerr2[i]/mm<<G4endl;
    //  G4cout <<hhh<<" "<<y[i]<<"  "<<ydif[i]/mm <<" "
    //	     <<yout1[i]/mm<<" "
    //       <<yout2[i]/mm<<" "<<yerr2[i]/mm<<G4endl;
}
///
///OneStepStepper
///
///
void OneStep(G4double y[],G4double dydx[],G4double hhh,
                 G4double yout1[],G4double yerr1[],
                  G4MagIntegratorStepper* pSt1)
{   

    pSt1->Stepper(y,dydx,hhh,yout1,yerr1);
    //G4cout << " One step with  stepper type 1"<<G4endl;
    //G4cout << " Step Length= "<<hhh/mm<<"mm"<<G4endl;
    //       G4cout << "Y    DyDx      Yout "<<G4endl;
    //G4cout.precision(10);
    //for (G4int i=0;i<6;i++)
    //G4cout << y[i]<<"  "<<dydx[i]<<"  "<<yout1[i]<<G4endl;

}
///
///OneStepStepper
///
 #ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
 #endif

int main(int argc,char** argv)
{
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;
 
  // set mandatory initialization classes
  //
 
   ExN01DetectorConstruction* detector = new ExN01DetectorConstruction;
   runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN01PhysicsList;
  runManager->SetUserInitialization(physics);
 
 #ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    //
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
 #endif

 // set mandatory user action class
  //
   G4VUserPrimaryGeneratorAction* gen_action = new ExN01PrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
 

 ExN01RunAction* run_action = new ExN01RunAction;
  runManager->SetUserAction(run_action);
  //
  ExN01EventAction* event_action = new ExN01EventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action =
    new ExN01SteppingAction(detector,event_action);
  runManager->SetUserAction(stepping_action);
  //
  //Mag Field
  //
   G4int type=4;
   G4ThreeVector Bconst(0.,0.,-2.*tesla);
   
   G4UniformMagField* pmagField;
   G4FieldManager* pfieldMgr;
   G4MagIntegratorStepper* pStepper1;
   G4MagIntegratorStepper* pStepper2;
   G4ChordFinder* pChordFinder;
   G4Mag_UsualEqRhs* pEquation;
   
   pmagField = new G4UniformMagField(Bconst);
  
   //// Get the global field manager 
   pfieldMgr= G4TransportationManager::GetTransportationManager()->
      GetFieldManager();
   // // Set this field to the global field manager 
   pfieldMgr->SetDetectorField(pmagField );
   // //stepper choise
    pEquation=new G4Mag_UsualEqRhs(pmagField);
    pStepper1=SetupStepper(type,pEquation);
    pChordFinder = new G4ChordFinder(pmagField,
   				       1.0e-2 * mm,
                                      pStepper1);
    pfieldMgr->SetChordFinder( pChordFinder );
      
   //
   //Propagator
   //
     G4PropagatorInField* pPropagator;
    pPropagator = G4TransportationManager::GetTransportationManager()->
       GetPropagatorInField ();
       pPropagator -> SetMinimumEpsilonStep( 1.0e-5*mm ) ;
    pPropagator -> SetMaximumEpsilonStep( 1.0e-5*mm ) ;
    G4cout << " Using values for "
           << " Min Eps = "  <<   pPropagator->GetMinimumEpsilonStep()
           << " and "
           << " MaxEps = " <<  pPropagator->GetMaximumEpsilonStep()
           << G4endl;
    pPropagator ->SetLargestAcceptableStep(5.0*cm);
   
   //
   //End MagField
 
   //

     // Initialize G4 kernel
  //
  runManager->Initialize();
// Get the pointer to the User Interface manager
  //
   G4UImanager * UI = G4UImanager::GetUIpointer();  
 
   if(argc==1)  // Define (G)UI terminal for interactive mode
   { 
     // G4UIterminal is a (dumb) terminal
     //
     G4UIsession * session = 0;
 #ifdef G4UI_USE_TCSH
       session = new G4UIterminal(new G4UItcsh);      
 #else
       session = new G4UIterminal();
 #endif
         UI->ApplyCommand("/control/execute vis.mac");  
	  
       //  UI->ApplyCommand("/gun/energy 1 GeV");   
	G4int numberOfEvent = 1;
	runManager->BeamOn(numberOfEvent);
        UI->ApplyCommand("/run/verbose 1"); 
       
       session->SessionStart();
     
     delete session;
   }
  else   // Batch mode
   { 
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
   }
 
   // Free the store: user actions, physics_list and detector_description are
   //                 owned and deleted by the run manager, so they should not
   //                 be deleted in the main() program !
 
  // Get the pointer to the UI manager and set verbosities
  //
  //G4UImanager* UI = G4UImanager::GetUIpointer();
  //UI->ApplyCommand("/run/verbose 1");
  //UI->ApplyCommand("/event/verbose 1");
  // UI->ApplyCommand("/tracking/verbose 1");

  // Start a run
  //
  //G4int numberOfEvent = 3;
  //runManager->BeamOn(numberOfEvent);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
#ifdef G4VIS_USE
   delete visManager;
 #endif
  delete runManager;

  return 0;
}


