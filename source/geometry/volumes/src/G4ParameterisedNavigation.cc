// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParameterisedNavigation.cc,v 1.1 1999-01-07 16:08:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4ParameterisedNavigation Implementation
//

#include "G4ParameterisedNavigation.hh"
G4ParameterisedNavigation::~G4ParameterisedNavigation()
{
#ifdef G4DEBUG_NAVIGATION
    cout << "G4ParameterisedNavigation::~G4ParameterisedNavigation() called."
	 << endl;
#endif
}

G4double G4ParameterisedNavigation::ComputeStep(const G4ThreeVector &localPoint,
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

	G4VPhysicalVolume *motherPhysical,*samplePhysical;
	G4VPVParameterisation *sampleParam;
	G4LogicalVolume *motherLogical;
	G4VSolid *motherSolid,*sampleSolid;
	G4ThreeVector sampleDirection;
	G4double ourStep=currentProposedStepLength,motherSafety,ourSafety;
	G4int sampleNo,blockedExitedReplicaNo=-1;

	G4bool initialNode,noStep;
	G4SmartVoxelNode *curVoxelNode;
	G4int curNoVolumes,contentNo;
	G4double voxelSafety;

// Replication data
	EAxis axis;
	G4int nReplicas;
	G4double width,offset;
	G4bool consuming;

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

	initialNode=true;
	noStep=true;

// By definition, parameterised volumes exist as first
// daughter of mother volume
	samplePhysical=motherLogical->GetDaughter(0);
	samplePhysical->GetReplicationData(axis,nReplicas,width,offset,consuming);
	fBList.Enlarge(nReplicas);
	fBList.Reset();


// Exiting normal optimisation
	if (exiting && (*pBlockedPhysical==samplePhysical) && validExitNormal)
		{
		if (localDirection.dot(exitNormal)>=kMinExitingNormalCosine)
			{
			assert(  (0 <= blockedReplicaNo)
			       &&(blockedReplicaNo<nReplicas));
// Block exited daughter replica; Must be on boundary => zero safety
			fBList.BlockVolume(blockedReplicaNo);
			ourSafety=0;
			}
		}

	exiting=false;
	entering=false;

	// sampleSolid=samplePhysical ->GetLogicalVolume() ->GetSolid();

	sampleParam=samplePhysical->GetParameterisation();

	do {

	curVoxelNode=fVoxelNode;
	curNoVolumes=curVoxelNode->GetNoContained();

	
	for (contentNo=curNoVolumes-1;contentNo>=0;contentNo--)
		{
		sampleNo=curVoxelNode->GetVolume(contentNo);
		if (!fBList.IsBlocked(sampleNo))
			{
			fBList.BlockVolume(sampleNo);
	                sampleSolid=sampleParam->ComputeSolid(sampleNo,
						 samplePhysical);
			sampleSolid->ComputeDimensions(sampleParam,
                                       	         sampleNo,
                                               	 samplePhysical);
			sampleParam->ComputeTransformation(sampleNo,
						 samplePhysical);
			samplePhysical->Setup(motherPhysical);
			G4AffineTransform sampleTf(samplePhysical->GetRotation(),
						   samplePhysical->GetTranslation());
			sampleTf.Invert();
			const G4ThreeVector samplePoint=sampleTf.TransformPoint(localPoint);
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
					blockedReplicaNo=sampleNo;
					}
	
				}
				    
		
			}
		}

	if (initialNode)
		{
		initialNode=false;
		voxelSafety=ComputeVoxelSafety(localPoint);
		if (voxelSafety<ourSafety)
			{
			ourSafety=voxelSafety;
			}

		if (currentProposedStepLength<ourSafety)
			{
//
// Guaranteed physics limited
//			
			noStep=false;
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
		}

	if (noStep)
		{
		noStep=LocateNextVoxel(localPoint,
				       localDirection,
				       ourStep);
		}

	} while (noStep);

	return ourStep;
}

G4double G4ParameterisedNavigation::ComputeSafety(const G4ThreeVector &localPoint,
			          const G4NavigationHistory &history,
		                  const G4double pProposedMaxLength )
{

	G4VPhysicalVolume *motherPhysical,*samplePhysical;
	G4VPVParameterisation *sampleParam;
	G4LogicalVolume *motherLogical;
	G4VSolid *motherSolid,*sampleSolid;
	G4double motherSafety,ourSafety;
	G4int sampleNo, curVoxelNodeNo;
	
	G4SmartVoxelNode *curVoxelNode;
	G4int curNoVolumes,contentNo;
	G4double voxelSafety;

// Replication data
	EAxis axis;
	G4int nReplicas;
	G4double width,offset;
	G4bool consuming;

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

// By definition, parameterised volumes exist as first
// daughter of mother volume
	samplePhysical=motherLogical->GetDaughter(0);
	samplePhysical->GetReplicationData(axis,nReplicas,width,offset,consuming);
	sampleParam=samplePhysical->GetParameterisation();

	// Calculate new VoxelNode of current point
        curVoxelNodeNo= G4int ( 
                          (localPoint(fVoxelAxis) -fVoxelHeader->GetMinExtent())
	                  / fVoxelSliceWidth
                        );
        curVoxelNode = fVoxelHeader->GetSlice(curVoxelNodeNo)->GetNode();

	curNoVolumes=curVoxelNode->GetNoContained();

	for (contentNo=curNoVolumes-1;contentNo>=0;contentNo--)
	    {
		sampleNo=curVoxelNode->GetVolume(contentNo);
			
		sampleSolid=sampleParam->ComputeSolid(sampleNo,
					 samplePhysical);
		sampleSolid->ComputeDimensions(sampleParam,
					 sampleNo,
					 samplePhysical);
		sampleParam->ComputeTransformation(sampleNo,
					 samplePhysical);
		G4AffineTransform sampleTf(samplePhysical->GetRotation(),
					   samplePhysical->GetTranslation());
		sampleTf.Invert();
		const G4ThreeVector samplePoint=sampleTf.TransformPoint(localPoint);
		const G4double sampleSafety=sampleSolid
				->DistanceToIn(samplePoint);
		if (sampleSafety<ourSafety)
			{
			ourSafety=sampleSafety;
			}
	    }

	// These must be current for ComputeVoxelSafety
        fVoxelNodeNo= curVoxelNodeNo;
        fVoxelNode =  curVoxelNode;

	voxelSafety=ComputeVoxelSafety(localPoint);
	if (voxelSafety<ourSafety)
	    {
		ourSafety=voxelSafety;
	    }

	return ourSafety;
}



















