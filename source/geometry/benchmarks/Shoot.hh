// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Shoot.hh,v 1.4 2000-12-12 08:23:08 medernac Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef SHOOT_HH
#define SHOOT_HH

#include "G4Timer.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

void Shoot(const G4int numShoot,
	   G4VPhysicalVolume *pTopNode,
	   const G4ThreeVector& pSource,
	   const G4ThreeVector& pVec)
{
  const G4double physStep=kInfinity;
  G4int i;
  G4double safety;
  G4double Step;
  G4Navigator myNav;
  G4Timer timer;
  G4ThreeVector partLoc;
  G4VPhysicalVolume *located;

  myNav.SetWorldVolume(pTopNode);

  timer.Start();

  for (i=numShoot;i>0;i--)
    {
      //      G4cout << "#Loop  " << i << G4endl ;
      
      partLoc=pSource;
      //located=myNav.LocateGlobalPointAndSetup(partLoc,false);
      located=myNav.LocateGlobalPointAndSetup(partLoc);
      while (located)
	{
	  /*
	  G4cout << "Loc = " << partLoc << " Vec = " << pVec << G4endl ;
	  G4cout << "Safety = " << safety << G4endl ;
	  */
	  Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
	  partLoc+=Step*pVec;
	  myNav.SetGeometricallyLimitedStep();
	  located=myNav.LocateGlobalPointAndSetup(partLoc);
	};
    }
  timer.Stop();
  //  G4cout << "Shots = " << numShoot << " " << timer << G4endl;
}

#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "PhysicalConstants.h"

void MagneticShoot(const G4int numShoot,
		   G4VPhysicalVolume *pTopNode,
		   const G4ThreeVector& pSource,
		   const G4ThreeVector& pVec,
		   const G4double fieldValue, // ** already in tesla **
		   const G4double DeltaChord) // ** already in mm **
{
  /** Setting up the Magnetic field **/

  G4UniformMagField magField (0.,0.,fieldValue);

  G4Navigator *myNav = G4TransportationManager::
    GetTransportationManager()-> GetNavigatorForTracking();
  myNav->SetWorldVolume(pTopNode);


  
  /* Field Properties */

  G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(&magField);
  
  /* Choose your stepper here */
  /* G4ClassicalRK4 is the default one */
  G4MagIntegratorStepper*  pStepper = new G4ClassicalRK4( fEquation );

  /*
    pStepper = new G4ExplicitEuler( fEquation );
    pStepper = new G4ImplicitEuler( fEquation );
    pStepper = new G4SimpleRunge( fEquation );
    pStepper = new G4SimpleHeum( fEquation );
    pStepper = new G4ClassicalRK4( fEquation );
    pStepper = new G4HelixExplicitEuler( fEquation );
    pStepper = new G4HelixImplicitEuler( fEquation );
    pStepper = new G4HelixSimpleRunge( fEquation );
    pStepper = new G4RKG3_Stepper( fEquation );
  */

  G4FieldManager* pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  
  pFieldMgr->SetDetectorField( &magField );

  //  G4ChordFinder* pChordFinder = new G4ChordFinder( &magField);
  G4ChordFinder* pChordFinder = new G4ChordFinder( &magField,DeltaChord,pStepper);

  pFieldMgr->SetChordFinder( pChordFinder );

  G4Timer timer;
  timer.Start();

  /*    
	G4PropagatorInField *pMagFieldPropagator = new G4PropagatorInField(myNav,pFieldMgr);
  */
  G4PropagatorInField *pMagFieldPropagator= G4TransportationManager::
    GetTransportationManager()-> GetPropagatorInField ();

  pChordFinder->SetChargeMomentumMass(  +1.,                    // charge in e+ units
					0.05 * proton_mass_c2,    // Momentum in Mev/c ?
					proton_mass_c2 );
  

  
  for (G4int i=numShoot;i>0;i--)
    {
      G4VPhysicalVolume *located;
      G4ThreeVector Vec = pVec ;
      G4ThreeVector partLoc = pSource ;	    
      /*
	G4cout << "#Loop  " << i << G4endl ;
	G4cout << "Loc = " << partLoc << " Vec = " << Vec << G4endl << G4endl ;
      */
      pFieldMgr->GetChordFinder()
	->SetChargeMomentumMass(// charge in e+ units
				+1,                  
				// Momentum in Mev/c 
				(0.5+i*0.1) * proton_mass_c2, 
				// Mass
				proton_mass_c2);

      
      
      located=myNav->LocateGlobalPointAndSetup(partLoc);
      
      while (located)
	{
	  const G4double physStep= kInfinity; // 2.5*mm ;
	  G4double safety = 1.0*m;
	  G4double Step = 0.0*m;
	  /*
	  G4cout << "Loc = " << partLoc << " Vec = " << Vec << G4endl ;
	  G4cout << "Safety = " << safety << G4endl ;
	  */
	  
	  Step=pMagFieldPropagator->ComputeStep(partLoc,Vec,physStep,safety);

	  myNav->SetGeometricallyLimitedStep();
	  
	  partLoc = pMagFieldPropagator->EndPosition();
	  Vec = pMagFieldPropagator->EndMomentumDir(); 
	  
	  located=myNav->LocateGlobalPointAndSetup(partLoc);
	};
    }
  timer.Stop();
  //  G4cout << "Shots = " << numShoot << " " << timer << G4endl;
}


void ShootVerbose(G4VPhysicalVolume *pTopNode,
		  const G4ThreeVector& pSource,
		  const G4ThreeVector& pVec)
{
  const G4double physStep=kInfinity;
  G4double safety,Step;
  G4Navigator myNav;
  G4ThreeVector partLoc;
  G4VPhysicalVolume *located=0;

  myNav.SetWorldVolume(pTopNode);

  partLoc=pSource;
  located=myNav.LocateGlobalPointAndSetup(partLoc);
  while (located)
    {
      Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
      G4cout << "Physical Location=" << located->GetName()
	     << " #" << located->GetCopyNo() << G4endl
	     << "   Step=" << Step << "  Safety=" << safety
	     << "  ---->" << G4endl;

      partLoc+=Step*pVec;
      myNav.SetGeometricallyLimitedStep();
      located=myNav.LocateGlobalPointAndSetup(partLoc);
    };
}

#endif












