// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Shoot.hh,v 1.1 1999-01-08 16:31:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef SHOOT_HH
#define SHOOT_HH

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4Timer.hh"
#include "G4ios.hh"

void Shoot(const G4int numShoot,
	   G4VPhysicalVolume *pTopNode,
	   const G4ThreeVector& pSource,
	   const G4ThreeVector& pVec)
{
    const G4double physStep=kInfinity;
    G4int i;
    G4double safety,Step;
    G4Navigator myNav;
    G4Timer timer;
    G4ThreeVector partLoc;
    G4VPhysicalVolume *located;

    myNav.SetWorldVolume(pTopNode);

    timer.Start();

    for (i=numShoot;i>0;i--)
	{
	    
	    partLoc=pSource;
	    located=myNav.LocateGlobalPointAndSetup(partLoc,false);
	    while (located)
		{
		    Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
		    partLoc+=Step*pVec;
		    myNav.SetGeometricallyLimitedStep();
		    located=myNav.LocateGlobalPointAndSetup(partLoc);
		};
	}
    timer.Stop();
    G4cout << "Shots = " << numShoot << " " << timer << endl;
}

void ShootVerbose(G4VPhysicalVolume *pTopNode,
		  const G4ThreeVector& pSource,
		  const G4ThreeVector& pVec)
{
    const G4double physStep=kInfinity;
    G4double safety,Step;
    G4Navigator myNav;
    G4ThreeVector partLoc;
    G4VPhysicalVolume *located;

    myNav.SetWorldVolume(pTopNode);

    partLoc=pSource;
    //    located=myNav.LocateGlobalPointAndSetup(partLoc,false);
    located=myNav.LocateGlobalPointAndSetup(partLoc,true);
    while (located)
	{
	    Step=myNav.ComputeStep(partLoc,pVec,physStep,safety);
	    G4cout << "Physical Location=" << located->GetName()
		 << " #" << located->GetCopyNo() << endl
	         << "   Step=" << Step << "  Safety=" << safety
		 << "  ---->" << endl;

	    partLoc+=Step*pVec;
	    myNav.SetGeometricallyLimitedStep();
	    located=myNav.LocateGlobalPointAndSetup(partLoc);
	};

}

#endif
