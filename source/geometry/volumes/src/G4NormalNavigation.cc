// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NormalNavigation.cc,v 1.3 2000-11-20 19:05:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4NormalNavigation Implementation
//

#include "G4NormalNavigation.hh"

G4double G4NormalNavigation::ComputeStep(const G4ThreeVector &localPoint,
			     const G4ThreeVector &localDirection,
			     const G4double currentProposedStepLength,
			     G4double &newSafety,
			     G4NavigationHistory &history,
			     G4bool &validExitNormal,
			     G4ThreeVector &exitNormal,
			     G4bool &exiting,
			     G4bool &entering,
                             G4VPhysicalVolume *(*pBlockedPhysical),
                             G4int &blockedReplicaNo)
{

	G4VPhysicalVolume *motherPhysical,*samplePhysical,*blockedExitedVol=0;
	G4LogicalVolume *motherLogical;
	G4VSolid *motherSolid;
	G4ThreeVector sampleDirection;
	G4double ourStep=currentProposedStepLength,motherSafety,ourSafety;
	G4int localNoDaughters,sampleNo;

	motherPhysical=history.GetTopVolume();
	motherLogical=motherPhysical->GetLogicalVolume();

	motherSolid=motherLogical->GetSolid();

//
// Compute mother safety
//
	motherSafety=motherSolid->DistanceToOut(localPoint);
	ourSafety=motherSafety; // Working isotropic safety
	
//
// Compute daughter safeties & intersections
//

// Exiting normal optimisation
	if (exiting&&validExitNormal)
		{
		if (localDirection.dot(exitNormal)>=kMinExitingNormalCosine)
			{
// Block exited daughter volume
			blockedExitedVol=*pBlockedPhysical;
			ourSafety=0;
			}
		}

	exiting=false;
	entering=false;

	localNoDaughters=motherLogical->GetNoDaughters();
	for (sampleNo=localNoDaughters-1;sampleNo>=0;sampleNo--)
		{
		samplePhysical=motherLogical->GetDaughter(sampleNo);
		if (samplePhysical!=blockedExitedVol)
			{
			samplePhysical->Setup(motherPhysical);
			G4AffineTransform sampleTf(samplePhysical->GetRotation(),
						   samplePhysical->GetTranslation());
			sampleTf.Invert();
			const G4ThreeVector samplePoint=sampleTf.TransformPoint(localPoint);
			const G4VSolid *sampleSolid=samplePhysical
						->GetLogicalVolume()
						->GetSolid();
			const G4double sampleSafety=sampleSolid
						->DistanceToIn(samplePoint);
			if (sampleSafety<ourSafety)
				{
			    	ourSafety=sampleSafety;
				}
			if (sampleSafety<=ourStep)
				{
				  sampleDirection=sampleTf.TransformAxis(localDirection);
				  const G4double sampleStep=sampleSolid
				               	->DistanceToIn(samplePoint,
							sampleDirection);
			    	if (sampleStep<=ourStep)
					{
					ourStep=sampleStep;
					entering=true;
					exiting=false;
					*pBlockedPhysical=samplePhysical;
					blockedReplicaNo=-1;
					}
				}
			}
		}
	if (currentProposedStepLength<ourSafety)
		{
//
// Guaranteed physics limited
//
		entering=false;
		exiting=false;
		*pBlockedPhysical=0;
		ourStep=kInfinity;
		}
	else
		{
//
// Compute mother intersection if required
//
		if (motherSafety<=ourStep)
			{
			G4double motherStep=motherSolid
					->DistanceToOut(localPoint,
							localDirection,
							true,
							&validExitNormal,
							&exitNormal);
			if (motherStep<=ourStep)
				{
				ourStep=motherStep;
				exiting=true;
				entering=false;
				if (validExitNormal)
					{
					const G4RotationMatrix *rot=motherPhysical->GetRotation();
					if (rot)
						{
						exitNormal*=rot->inverse();
						}
					}
				}
			else
				{
				validExitNormal=false;
				}
			}
		}
	newSafety=ourSafety;
	return ourStep;
}

G4double G4NormalNavigation::ComputeSafety(const G4ThreeVector &localPoint,
			     const G4NavigationHistory &history,
			     const G4double)
{
     G4VPhysicalVolume *motherPhysical,*samplePhysical;
     G4LogicalVolume *motherLogical;
     G4VSolid *motherSolid;
     G4double motherSafety,ourSafety;
     G4int localNoDaughters,sampleNo;

     motherPhysical=history.GetTopVolume();
     motherLogical=motherPhysical->GetLogicalVolume();
     motherSolid=motherLogical->GetSolid();

   //
   // Compute mother safety
   //
     motherSafety=motherSolid->DistanceToOut(localPoint);
     ourSafety=motherSafety; // Working isotropic safety

   //
   // Compute daughter safeties 
   //
     localNoDaughters=motherLogical->GetNoDaughters();
     for (sampleNo=localNoDaughters-1;sampleNo>=0;sampleNo--)
	{
	    samplePhysical=motherLogical->GetDaughter(sampleNo);

	    samplePhysical->Setup(motherPhysical);
	    G4AffineTransform sampleTf(samplePhysical->GetRotation(),
				       samplePhysical->GetTranslation());
	    sampleTf.Invert();
	    const G4ThreeVector samplePoint=sampleTf.TransformPoint(localPoint);
	    const G4VSolid *sampleSolid=samplePhysical
				    ->GetLogicalVolume()
				    ->GetSolid();
	    const G4double sampleSafety=sampleSolid
				    ->DistanceToIn(samplePoint);
	    if (sampleSafety<ourSafety)
		{
		    ourSafety=sampleSafety;
		}
	}

     return ourSafety;
}
