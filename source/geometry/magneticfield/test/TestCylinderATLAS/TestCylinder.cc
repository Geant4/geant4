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
// $Id: exampleN01.cc,v 1.6 2006/06/29 17:47:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01-ref $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
// --------------------------------------------------------------

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4HelixMixedStepper.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4Timer.hh"
#include "G4ios.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4GeometryManager.hh" 
#include "G4GeometryTolerance.hh" 
#include <fstream>
////Root for equation of second order
void root2(G4double a,G4double b, G4double c,G4double* x1,G4double* x2)
{
  G4double D;
  if(a==0.)G4cout<<"Root exeption : a==0"<<G4endl;
  D=b*b-4.*a*c;
  G4cout<<"Deter="<<D<<" a="<<a<<" b="<<b<<" c="<<c<<G4endl;
  if(D<0.00000001*mm){
    G4cout<<"No solution in Root:D<0"<<G4endl;
    *x1=0.;*x2=0.;
  }
  else {
    *x1=(-b-std::sqrt(D))/(2.*a);
    *x2=(-b+std::sqrt(D))/(2.*a);

  }
  G4cout<<"in root2 solution="<<*x1<<"  "<<*x2<<G4endl;
}
////End Root
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
      case 11: pStepper = new G4HelixMixedStepper( pE );    break;
      default: pStepper = 0;
    }
  return pStepper;
}
///
///EndSetupStepper
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
  // G4int type=9;
      G4int type=4;
   G4ThreeVector Bconst(0.,0.,-2.*tesla);
   
   G4UniformMagField* pmagField;
   G4FieldManager* pfieldMgr;
   G4MagIntegratorStepper* pStepper1;
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
    pfieldMgr->SetDeltaIntersection(1e-3*mm);
    //for many intersection  pfieldMgr->SetDeltaIntersection(1e-3*mm);
     G4cout<<"Delta= "<<pfieldMgr->GetDeltaIntersection()<<G4endl; 
    //
   //Propagator
   //
     G4PropagatorInField* pPropagator;
    pPropagator = G4TransportationManager::GetTransportationManager()->
       GetPropagatorInField ();
       pPropagator -> SetMinimumEpsilonStep( 1.0e-5 ) ;
    pPropagator -> SetMaximumEpsilonStep( 1.0e-3 ) ;
    G4cout << " Using values for "
           << " Min Eps = "  <<   pPropagator->GetMinimumEpsilonStep()
           << " and "
           << " MaxEps = " <<  pPropagator->GetMaximumEpsilonStep()
           << G4endl;
    // pPropagator ->SetLargestAcceptableStep(10.*mm);
    // pPropagator ->SetLargestAcceptableStep(3.0*mm);
   
   //
   //End MagField
  
    G4double R1=70.136*cm;
    G4double R2=70.1381368117*cm;
    G4double x0=1*mm;
    G4double y0=0.;
    G4double a,b,aa,bb,cc;
    G4double x1,y1,x2,y2;
    if(y0==0.){
      x1=(R1*R1+x0*x0-R2*R2)/(2.*x0);
      x2=x1;
      if((R2*R2-(x1-x0)*(x1-x0))<0.){
	y1=0.;y2=0;
      }
      else{
      y1=-std::sqrt(R2*R2-(x1-x0)*(x1-x0));
      y2=std::sqrt(R2*R2-(x1-x0)*(x1-x0));
      }
    }
    else{
      a=(R1*R1-R2*R2+x0*x0+y0*y0)/(2.*y0);
      b=x0/y0;
      aa=1+b*b;;
      bb=-2.*a*b;
      cc=a*a-R1*R1;
      root2(aa,bb,cc,&x1,&x2);
      y1=a-b*x1;y2=a-b*x2;
    }
    G4cout.precision(15);
    G4cout<<"Intersection Points :"<<G4endl;
    G4cout<<"(x1,y1)="<<x1<<" , "<<y1<<G4endl;
    G4cout<<"(x2,y2)="<<x2<<" , "<<y2<<G4endl;
    //
    //End Exact Intersection Calculation

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

#ifdef G4VIS_USE
   delete visManager;
 #endif
  delete runManager;

  return 0;
}


