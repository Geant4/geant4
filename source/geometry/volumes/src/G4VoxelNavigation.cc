// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VoxelNavigation.cc,v 1.2 1999-01-08 11:23:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4VoxelNavigation Implementation
//
// $ Id: $
//
// Modified by:
//   J. Apostolakis,  29 Apr 98 Fixed error in LocateNextVoxel that
//                              ignored voxels at lower levels 

#include "G4VoxelNavigation.hh"
#include <iomanip.h>

#define EBUG_VOXEL_NAV 1

G4double G4VoxelNavigation::ComputeStep(const G4ThreeVector &localPoint,
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

	G4bool initialNode,noStep;
	G4SmartVoxelNode *curVoxelNode;
	G4int curNoVolumes,contentNo;
	G4double voxelSafety;

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

	fBList.Enlarge(localNoDaughters);
	fBList.Reset();

	initialNode=true;
	noStep=true;
        G4int countStep= 0;

	do {

	curVoxelNode=fVoxelNode;
	curNoVolumes=curVoxelNode->GetNoContained();

	
	for (contentNo=curNoVolumes-1;contentNo>=0;contentNo--)
		{
		sampleNo=curVoxelNode->GetVolume(contentNo);
		if (!fBList.IsBlocked(sampleNo))
			{
			fBList.BlockVolume(sampleNo);
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
#ifdef EBUG_VOXEL_NAV
                if (countStep == 0)
		  {
		    G4cout << "**** G4VoxelNavigation::ComputeStep ****" << endl;
		    G4cout << " Location: "  << localPoint     << " " 
			   << " Direction: " << localDirection << " " << endl;
		    printHeaderLocNxVox();
		  }
		G4cout << setw(12) << "LocNxVox bef:" 
		       << setw(6) << countStep << " ";
		G4cout << setw(4) <<  "    " ; 
		printStateLocNxVox();
#endif

		noStep=LocateNextVoxel(localPoint,
				       localDirection,
				       ourStep);
#ifdef EBUG_VOXEL_NAV
		G4cout << setw(12) << "LocNxVox ret:" 
		       << setw(4)  << noStep 
		       << setw(6)  << countStep << " ";
		printStateLocNxVox();

                countStep++;
#endif
	  }

	} while (noStep);

	return ourStep;
}


// Compute safety from specified point to voxel boundaries
// using already located point
// o collected boundaries for most derived level
// o adjacent boundaries for previous levels

G4double G4VoxelNavigation::ComputeVoxelSafety(const G4ThreeVector&localPoint) const
	{
	G4SmartVoxelHeader *curHeader;
	G4double voxelSafety,curNodeWidth;
	G4double curNodeOffset,minCurCommonDelta,maxCurCommonDelta;
	G4int minCurNodeNoDelta,maxCurNodeNoDelta;
	G4int localVoxelDepth,curNodeNo;
	EAxis curHeaderAxis;

	localVoxelDepth=fVoxelDepth;

	curHeader=fVoxelHeaderStack(localVoxelDepth);
	curHeaderAxis=fVoxelAxisStack(localVoxelDepth);
	curNodeNo=fVoxelNodeNoStack(localVoxelDepth);
	curNodeWidth=fVoxelSliceWidthStack(localVoxelDepth);
	
// Compute linear intersection distance to boundaries of max/min
// to collected nodes at current level
	curNodeOffset=curNodeNo*curNodeWidth;
	maxCurNodeNoDelta=fVoxelNode->GetMaxEquivalentSliceNo()-curNodeNo;
	minCurNodeNoDelta=curNodeNo-fVoxelNode->GetMinEquivalentSliceNo();
	minCurCommonDelta=localPoint(curHeaderAxis)
	    -curHeader->GetMinExtent()
	    -curNodeOffset;
	maxCurCommonDelta=curNodeWidth-minCurCommonDelta;

	if (minCurNodeNoDelta<maxCurNodeNoDelta)
		{
		voxelSafety=minCurNodeNoDelta*curNodeWidth;
		voxelSafety+=minCurCommonDelta;
	    	}
	else if (maxCurNodeNoDelta<minCurNodeNoDelta)
	    	{
		voxelSafety=maxCurNodeNoDelta*curNodeWidth;
		voxelSafety+=maxCurCommonDelta;
	    	}
	else // (maxCurNodeNoDelta == minCurNodeNoDelta)
	    	{
		voxelSafety=minCurNodeNoDelta*curNodeWidth;
		voxelSafety+=min(minCurCommonDelta,maxCurCommonDelta);
	    	}

// Compute isotropic safety to boundaries of previous levels
// [NOT to collected boundaries]
	while (localVoxelDepth>0&&voxelSafety>0)
	    	{
		localVoxelDepth--;

		curHeader=fVoxelHeaderStack(localVoxelDepth);
		curHeaderAxis=fVoxelAxisStack(localVoxelDepth);
		curNodeNo=fVoxelNodeNoStack(localVoxelDepth);
		curNodeWidth=fVoxelSliceWidthStack(localVoxelDepth);
		curNodeOffset=curNodeNo*curNodeWidth;
		minCurCommonDelta=localPoint(curHeaderAxis)
		    -curHeader->GetMinExtent()
		    -curNodeOffset;
		maxCurCommonDelta=curNodeWidth-minCurCommonDelta;
		
		if (minCurCommonDelta<voxelSafety)
			{
			voxelSafety=minCurCommonDelta;
		    	}
		if (maxCurCommonDelta<voxelSafety)
		    	{
			voxelSafety=maxCurCommonDelta;
		    	}
	    	}

	if (voxelSafety<0)
		{
		voxelSafety=0;
		}

	return voxelSafety;
	}


// Find the next voxel from the current voxel and point in the specified
// direction
//
// Return false if all voxels considered
//              [current Step ends inside same voxel or leaves all voxels]
//        true  otherwise
//              [the information on the next voxel is put into the set of
//                fVoxel* variables & "stacks" ] 
//
// 
G4bool G4VoxelNavigation::LocateNextVoxel(const G4ThreeVector& localPoint,
					   const G4ThreeVector& localDirection,
					   const G4double currentStep)
{
	G4SmartVoxelHeader *workHeader,*newHeader;
	G4SmartVoxelProxy *newProxy;
	G4SmartVoxelNode *newVoxelNode;
	G4ThreeVector targetPoint,voxelPoint;
	G4double workNodeWidth,workMinExtent,workCoord;
	G4double minVal,maxVal,newDistance;
	G4double newHeaderMin,newHeaderNodeWidth;
	G4int depth, newDepth,workNodeNo,newNodeNo,newHeaderNoSlices;
	EAxis workHeaderAxis,newHeaderAxis;
	G4bool isNewVoxel=false;
	
	G4double currentDistance= currentStep;

#ifdef EBUG_VOXEL_NAV_MORE
        G4cout << "*** Entered G4VoxelNavigation::LocateNextVoxel ***" << endl;
	G4cout << " Arguments localPoint, localDirection are " 
	       << localPoint << localDirection << endl; 
        G4cout << " Voxel Depth = " << fVoxelDepth << endl;
#endif 

// Determine if end of Step within current voxel
	for (depth=0;depth<fVoxelDepth;depth++)
	    {
		targetPoint=localPoint+localDirection*currentDistance;
	        newDistance= currentDistance;
		workHeader=fVoxelHeaderStack(depth);
		workHeaderAxis=fVoxelAxisStack(depth);
		workNodeNo=fVoxelNodeNoStack(depth);
		workNodeWidth=fVoxelSliceWidthStack(depth);
		workMinExtent=workHeader->GetMinExtent();

		// #define EBUG_NEXT_VOXEL 1

#ifdef EBUG_NEXT_VOXEL
		G4cout << "WorkHeader " << endl ;
		G4cout << "Depth Axis NodeNo  Width    Extent  WorkCoord     minVal      maxVal " << endl;
		G4cout << setw(5) << depth 
		       << setw(5) << workHeaderAxis  
		       << setw(5) << workNodeNo << "   "
		       << setw(8) << workNodeWidth
		       << setw(8) << workMinExtent << " " ;
#endif

		workCoord=targetPoint(workHeaderAxis);
		minVal=workMinExtent+workNodeNo*workNodeWidth;
#ifdef EBUG_NEXT_VOXEL
		G4cout << setw(8) << workCoord << " " 
		       << setw(8) << minVal;
#endif
		if (minVal<=workCoord+kCarTolerance*0.5)
		    {
			maxVal=minVal+workNodeWidth;
#ifdef EBUG_NEXT_VOXEL
			G4cout << setw(8) << maxVal << endl;
#endif
			if (maxVal<=workCoord-kCarTolerance*0.5)
			  {
#ifdef EBUG_NEXT_VOXEL
			    G4cout << "Considering next voxel in +direction" << endl;
#endif
			    // Ensure that it is moving in negative direction
			    // (not simply within kCarTolerance of the surface)
                            //  Enforces a logical decision on the surface
			    if( localDirection(workHeaderAxis) > 0 ) 
			      {
				newNodeNo=workNodeNo+1;
				newHeader=workHeader;
				newDistance=(maxVal-localPoint(workHeaderAxis))/localDirection(workHeaderAxis);
				isNewVoxel=true;
				newDepth= depth;
			      }
			  }
		    }
		else
		    {
#ifdef EBUG_NEXT_VOXEL
		      G4cout << endl;
		      G4cout << "Considering next voxel in -direction" << endl;
#endif
                      // Ensure that it is moving in negative direction
		      //  (not simply within kCarTolerance of the surface)
		      if( localDirection(workHeaderAxis) < 0 ) 
                        {
				newNodeNo=workNodeNo-1;
				newHeader=workHeader;
				newDistance=(minVal-localPoint(workHeaderAxis))/localDirection(workHeaderAxis);
				isNewVoxel=true;
				newDepth= depth;
			}
		    }

	        if( newDistance >= kCarTolerance * 0.5 ) 
		  {
		    currentDistance= newDistance;
		  }
                else if ( newDistance >= - kCarTolerance * 0.5 ) 
                  {
                    currentDistance = newDistance = 0; 
		  }
#ifdef EBUG_VOXEL_NAV
		G4cout << " Depth = " << depth; 
		G4cout << " New distance = " << newDistance 
		       << " Current dist = " << currentDistance; 
		G4cout << endl;
		// assert ( newDistance >= 0 );
		// if( newDistance < 0 )  
		//  G4Exception(" ERROR in G4VoxelNavigation::LocateNextVoxel  newDistance is negative.");
#endif
	    }

#ifdef CHECK_PRE_CONDITIONS
        minVal=workMinExtent+fVoxelNode->GetMinEquivalentSliceNo()*workNodeWidth;
        assert( (minVal-localPoint(workHeaderAxis)) * localDirection(workHeaderAxis) >= 0 );        
#endif

	targetPoint=localPoint+localDirection*currentDistance;
// Check if end of Step within collected boundaries of current voxel
	depth=fVoxelDepth;
	    {
		workHeader=fVoxelHeaderStack(depth);
		workHeaderAxis=fVoxelAxisStack(depth);
		workNodeNo=fVoxelNodeNoStack(depth);
		workNodeWidth=fVoxelSliceWidthStack(depth);
		workMinExtent=workHeader->GetMinExtent();
		    
		workCoord=targetPoint(workHeaderAxis);
		minVal=workMinExtent+fVoxelNode->GetMinEquivalentSliceNo()*workNodeWidth;
		if (minVal<=workCoord+kCarTolerance*0.5)
		    {
			maxVal=workMinExtent+(fVoxelNode->GetMaxEquivalentSliceNo()+1)*workNodeWidth;
			if (maxVal<=workCoord-kCarTolerance*0.5)
			    {
			    // Ensure that it is moving in negative direction
			    // (not simply within kCarTolerance of the surface)
                            //  Enforces a logical decision on the surface
			    if( localDirection(workHeaderAxis) > 0 ) 
			      {
				newNodeNo=fVoxelNode->GetMaxEquivalentSliceNo()+1;
				newHeader=workHeader;
				newDistance=(maxVal-localPoint(workHeaderAxis))/localDirection(workHeaderAxis);
				isNewVoxel=true;
				newDepth= depth;
			      }
			    }
		    }
		else
		    {
                      // Ensure that it is moving in negative direction
		      //  (not simply within kCarTolerance of the surface)
		      if( localDirection(workHeaderAxis) < 0 ) 
                        {
				newNodeNo=fVoxelNode->GetMinEquivalentSliceNo()-1;
				newHeader=workHeader;
				newDistance=(minVal-localPoint(workHeaderAxis))/localDirection(workHeaderAxis);
				isNewVoxel=true;
				newDepth= depth;
			}
		    }

	        currentDistance= newDistance;

#ifdef EBUG_VOXEL_NAV
		G4cout << " Depth = " << depth; 
		G4cout << " New distance = " << currentDistance ; 
		G4cout << "  -- full depth, ie equivalent nodes used" << endl;
		// assert ( newDistance >= 0 );
		if( newDistance < 0 )  
		  G4Exception(" ERROR in G4VoxelNavigation::LocateNextVoxel  newDistance is negative.");
#endif
	    }

	if (isNewVoxel)
	    {
// Compute new voxel & adjust voxel stack
//
// newNodeNo=Candidate node no at 
// newDepth =refinement depth of crossed voxel boundary
// newHeader=Header for crossed voxel
// newDistance=distance to crossed voxel boundary (along the track)
//
		if (newNodeNo<0||newNodeNo>=newHeader->GetNoSlices())
		    {
// Leaving mother volume
			isNewVoxel=false;
		    }
		else
		    {
// Compute intersection point on the least refined voxel boundary that is Hit
			voxelPoint=localPoint+localDirection*newDistance;
			fVoxelNodeNoStack(newDepth)=newNodeNo;
			fVoxelDepth=newDepth;
			newVoxelNode=0;

			while (!newVoxelNode)
			    {
				newProxy=newHeader->GetSlice(newNodeNo);
				if (newProxy->IsNode())
				    {
					newVoxelNode=newProxy->GetNode();
				    }
				else
				    {
					fVoxelDepth++;
				        newHeader=newProxy->GetHeader();
					newHeaderAxis=newHeader->GetAxis();
					newHeaderNoSlices=newHeader->GetNoSlices();
					newHeaderMin=newHeader->GetMinExtent();
					newHeaderNodeWidth=(newHeader->GetMaxExtent()-newHeaderMin)/newHeaderNoSlices;
					newNodeNo=G4int ((voxelPoint(newHeaderAxis)-newHeaderMin)/newHeaderNodeWidth);
// Rounding protection
					if (newNodeNo<0)
					    {
						newNodeNo=0;
					    }
					else if (newNodeNo>=newHeaderNoSlices)
					    {
						newNodeNo=newHeaderNoSlices-1;
					    }
// Stack info for stepping
					fVoxelAxisStack(fVoxelDepth)=newHeaderAxis;
					fVoxelNoSlicesStack(fVoxelDepth)=newHeaderNoSlices;
					fVoxelSliceWidthStack(fVoxelDepth)=newHeaderNodeWidth;
					fVoxelNodeNoStack(fVoxelDepth)=newNodeNo;
					fVoxelHeaderStack(fVoxelDepth)=newHeader;
				    }
		
			    }
			fVoxelNode=newVoxelNode;
		    }
	    }

#ifdef EBUG_NEXT_VOXEL
        G4cout << "*** Exiting G4VoxelNavigation::LocateNextVoxel ***" << endl;
#endif
	return isNewVoxel;				
}

//-----------------------------------------------------------------------------


//  Calculate the isotropic distance to the nearest boundary from the
// specified point in the local coordinate system. 
// The localpoint utilised must be within the current volume.

G4double G4VoxelNavigation::ComputeSafety(const G4ThreeVector &localPoint,
			       const G4NavigationHistory &history,
    		               const G4double pMaxLength )
{

     G4VPhysicalVolume *motherPhysical,*samplePhysical;
     G4LogicalVolume *motherLogical;
     G4VSolid *motherSolid;
     G4double motherSafety,ourSafety;
     G4int localNoDaughters,sampleNo;

     G4SmartVoxelNode *curVoxelNode;
     G4int curNoVolumes,contentNo;
     G4double voxelSafety;

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

//
//  Look only inside the current Voxel only (in the first version).
//

     curVoxelNode=fVoxelNode;
     curNoVolumes=curVoxelNode->GetNoContained();

     for (contentNo=curNoVolumes-1;contentNo>=0;contentNo--)
	{
	   sampleNo=curVoxelNode->GetVolume(contentNo);
	   samplePhysical=motherLogical->GetDaughter(sampleNo);

	   samplePhysical->Setup(motherPhysical);
	   G4AffineTransform sampleTf(samplePhysical->GetRotation(),
				      samplePhysical->GetTranslation());
	   sampleTf.Invert();
	   const G4ThreeVector samplePoint=sampleTf.TransformPoint(localPoint);
	   const G4VSolid *sampleSolid=samplePhysical ->GetLogicalVolume()
						      ->GetSolid();
	   const G4double sampleSafety=sampleSolid
					->DistanceToIn(samplePoint);
	   if (sampleSafety<ourSafety)
		{
		   ourSafety=sampleSafety;
		}		
		  
	}

    voxelSafety=ComputeVoxelSafety(localPoint);
    if (voxelSafety<ourSafety)
	{
	    ourSafety=voxelSafety;
	}
    
    return ourSafety;
}

G4VoxelNavigation::~G4VoxelNavigation()
{
#ifdef G4DEBUG_NAVIGATION
  cout << "G4VoxelNavigation::~G4VoxelNavigation() called." << endl;
#endif
}

void G4VoxelNavigation::printHeaderLocNxVox()
{
  G4cout << setw(12) << "-Print From-" 
	 << setw(6)  << "Count "       << " "  
	 << setw(4)  << "bool"         
	 << setw(12) << "Voxel:Depth"  << " " 
	 << setw(10) << "NodeNumbers"  << " " 
	 << setw( 6) << "AxisNo" << " " 
	 << endl;
  //  G4cout << " Before calling LocateNextVoxel: " 
  //     << endl;
}

void G4VoxelNavigation::printStateLocNxVox()
{
  G4cout <<  setw(10) << fVoxelDepth << " ";
  G4int depth;
  for (depth=0;depth<=fVoxelDepth;depth++)
     G4cout <<  setw(3) << fVoxelNodeNoStack(depth) << " "; 
  for (depth=fVoxelDepth;depth<=3;depth++)
     G4cout <<  setw(3) << "    " << " " ; 
  G4cout << "    " ; 
  for (depth=0;depth<=fVoxelDepth;depth++)
     G4cout <<  setw(1) << fVoxelAxisStack(depth) << " "; 
  for (depth=fVoxelDepth;depth<=3;depth++)
     G4cout <<  setw(1) << " " << " " ; 
  G4cout << " " ; 
  G4cout << endl;
}
