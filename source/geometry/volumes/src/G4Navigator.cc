// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Navigator.cc,v 1.3 1999-01-29 18:38:57 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Navigator Implementation  Paul Kent July 95/96

#include "G4Navigator.hh"
#include "G4ios.hh"
#include <iomanip.h>

G4Navigator::G4Navigator() :
  fWasLimitedByGeometry(false),
  fTopPhysical(0),
  fVerbose(0)
{
  ResetStackAndState();
}
    
G4Navigator::~G4Navigator()
{;}

// Set the world (`topmost') volume
void G4Navigator::SetWorldVolume(G4VPhysicalVolume* pWorld)
{
// Setup the volume
  pWorld->Setup(0);	// No mother since world volume
  if (!(pWorld->GetTranslation()==G4ThreeVector(0,0,0)))
    {
      G4Exception ("G4Navigator::SetWorldVolume - Must be centred on origin");
    }
  const G4RotationMatrix* rm=pWorld->GetRotation();
  if (rm&&(!rm->isIdentity()))
    {
      G4Exception ("G4Navigator::SetWorldVolume - Must not be rotated");
    }
  fTopPhysical=pWorld;
  fHistory.SetFirstEntry(pWorld);
}


// define DEBUG_HIST 1

// Locate the point in the hierarchy return 0 if outside
//
//  ( The direction is required only if we are on an edge shared by
//     two or more surfaces. )
// 
G4VPhysicalVolume* 
G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector& globalPoint,
				       const G4ThreeVector* pGlobalDirection,
				       const G4bool relativeSearch)
{
  G4bool notKnownContained=true,noResult;
  G4VPhysicalVolume *targetPhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *targetSolid;
  G4ThreeVector localPoint;
  EInside insideCode;

#ifdef DEBUG_HIST
  G4cerr << "Upon entering LocateGlobalPointAndSetup " << endl;
  G4cerr << "  History = " << endl << fHistory << endl << endl;
#endif
#ifdef G4VERBOSE
  if( fVerbose > 0 ) 
    {
      cout << "*** G4Navigator::LocateGlobalPointAndSetup: ***" << endl; 
      cout.precision(8);
      cout << " I was called with the following arguments: " << endl
	   << " Globalpoint = " << globalPoint << endl
	   << " relativeSearch   = " <<  relativeSearch  << endl;
       //       << " = " << << endl
      cout << " Upon entering my state is: " << endl;
      PrintState();
    }
#endif

  if (!relativeSearch)
    {
      ResetStackAndState();
    }
  else
    {
      if (fWasLimitedByGeometry)
	{
	  fWasLimitedByGeometry=false;
	  fEnteredDaughter=fEntering;   // Remember
          fExitedMother= fExiting;      // Remember
          if (fExiting)
	    {
	      if (fHistory.GetDepth())
		{
		  fBlockedPhysicalVolume=fHistory.GetTopVolume();
		  fBlockedReplicaNo=fHistory.GetTopReplicaNo();
		  fHistory.BackLevel();
		}
	      else
		{
		  // Have exited world volume
		  return 0;
		}

	      // A fix for the case where a volume is "entered" at an edge
	      //   and a coincident surface exists outside it.
	      //  This stops it from exiting further volumes and cycling
	      if( fLastStepWasZero )
		{ 
		  fExiting= false;
		}
	    }
	  else if (fEntering)
	    {
	      G4VPhysicalVolume *curPhysical=fHistory.GetTopVolume();
	      switch (VolumeType(fBlockedPhysicalVolume))
		{
		case kNormal:
		  fBlockedPhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fBlockedPhysicalVolume);
		  break;
		case kReplica:
		  freplicaNav.ComputeTransformation(fBlockedReplicaNo,
						    fBlockedPhysicalVolume);
		  fBlockedPhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fBlockedPhysicalVolume,
				    kReplica,
				    fBlockedReplicaNo);
		  fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
		  
		  break;
		case kParameterised:
		  G4VSolid *pSolid;
		  // G4VSolid *pSolid=fBlockedPhysicalVolume->
		  // GetLogicalVolume()-> GetSolid();
		  G4VPVParameterisation *pParam=fBlockedPhysicalVolume->
		    GetParameterisation();
		  pSolid= pParam->ComputeSolid(fBlockedReplicaNo,
						fBlockedPhysicalVolume);
		  pSolid->ComputeDimensions(pParam,
					    fBlockedReplicaNo,
					    fBlockedPhysicalVolume);
		  pParam->ComputeTransformation(fBlockedReplicaNo,
						fBlockedPhysicalVolume);
		  fBlockedPhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fBlockedPhysicalVolume,
				    kParameterised,
				    fBlockedReplicaNo);
		  fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
		  // Set the correct solid and material in Logical Volume
		  G4LogicalVolume *pLogical;
		  pLogical= fBlockedPhysicalVolume->GetLogicalVolume();
		  pLogical->SetSolid( pSolid );
		  pLogical->SetMaterial( 
			     pParam->ComputeMaterial(fBlockedReplicaNo, 
						     fBlockedPhysicalVolume));
		  break;
		}
	      fEntering=false;
	      fBlockedPhysicalVolume=0;
	      localPoint=fHistory.GetTopTransform().TransformPoint(globalPoint);
	      notKnownContained=false;
	    }
	}
      else
	{
	  fBlockedPhysicalVolume=0;
	  fEntering=false;
	  fEnteredDaughter=false;  // Full Step was not taken, did not enter
	  fExiting=false;
          fExitedMother=false;     // Full Step was not taken, did not exit
	}
    }
//
// Search from top of history up through geometry until
// containing volume found:
//
// If on 
// o OUTSIDE - Back up level, not/no longer exiting volumes
// o SURFACE and EXITING - Back up level, setting new blocking no.s
// else
// o containing volume found
//
  while (notKnownContained)
    {
      if (fHistory.GetTopVolumeType()!=kReplica)
	{
	  targetSolid=fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid();
	  localPoint=fHistory.GetTopTransform().TransformPoint(globalPoint);
	  insideCode=targetSolid->Inside(localPoint);
	}
      else
	{
	  insideCode=freplicaNav.BackLocate(fHistory,globalPoint,localPoint,fExiting,notKnownContained);
	  // !CARE! if notKnownContained returns false then the point is within
	  // the containing placement volume of the replica(s). If insidecode
	  // will result in the history being backed up one level, then the
	  // local point returned is the point in the system of this new level
	}

      if (insideCode==kOutside)
	{
	  if (fHistory.GetDepth())
	    {
	      fBlockedPhysicalVolume=fHistory.GetTopVolume();
	      fBlockedReplicaNo=fHistory.GetTopReplicaNo();
	      fHistory.BackLevel();
	      fExiting=false;
	    }
	  else
	    {
	      // Have exited world volume
	      return 0;
	    }			
	}
      else if (insideCode==kSurface&&fExiting)
	{
	  if (fHistory.GetDepth())
	    {
	      fBlockedPhysicalVolume=fHistory.GetTopVolume();
	      fBlockedReplicaNo=fHistory.GetTopReplicaNo();
	      fHistory.BackLevel();
	      // Still on surface but exited volume not necessarily convex
	      fValidExitNormal=false;
	    }
	  else
	    {
	      // Have exited world volume
	      return 0;
	    }			
	}
      else
	{
	  notKnownContained=false;
	}
    }
  
//
// Search downwards until deepest containing volume found,
// blocking fBlockedPhysicalVolume/BlockedReplicaNum
//
// 3 Cases:
//
// o Parameterised daughters
//   =>Must be one G4PVParameterised daughter & voxels
// o Positioned daughters & voxels
// o Positioned daughters & no voxels

  noResult=true;  // noResult should be renamed to 
                  //  something like enteredLevel, as that is its meaning.
  do
    {
      
// Determine `type' of current mother volume
      targetPhysical=fHistory.GetTopVolume();
      targetLogical=targetPhysical->GetLogicalVolume();
      switch(CharacteriseDaughters(targetLogical))
	{
	case kNormal:
	  if (targetLogical->GetVoxelHeader())
	    {
	      noResult=fvoxelNav.LevelLocate(fHistory,
					     fBlockedPhysicalVolume,
					     fBlockedReplicaNo,
					     globalPoint,
					     pGlobalDirection,
					     fLocatedOnEdge,
					     localPoint);
	    }
	  else
	    {
	      noResult=fnormalNav.LevelLocate(fHistory,
					      fBlockedPhysicalVolume,
					      fBlockedReplicaNo,
					      globalPoint,
					      pGlobalDirection,
					      fLocatedOnEdge,
					      localPoint);
	    }
	  break;
	case kReplica:
	  noResult=freplicaNav.LevelLocate(fHistory,
					   fBlockedPhysicalVolume,
					   fBlockedReplicaNo,
					   globalPoint,
					   pGlobalDirection,
					   fLocatedOnEdge,
					   localPoint);
	  break;
	case kParameterised:
	  noResult=fparamNav.LevelLocate(fHistory,
					 fBlockedPhysicalVolume,
					 fBlockedReplicaNo,
					 globalPoint,
					 pGlobalDirection,
					 fLocatedOnEdge,
					 localPoint);
	  break;
	}

      // LevelLocate returns true if it finds a daughter volume 
      //  in which globalPoint is inside (or on the surface).

      if (noResult)
	{
	  // The blocked volume is no longer valid - it was for another level
	  fBlockedPhysicalVolume= 0;
          fBlockedReplicaNo= -1;
	}
    } while (noResult);

  fLastLocatedPointLocal=localPoint;

#ifdef G4VERBOSE
  if( fVerbose > 0 )  PrintState();
    
  if( fVerbose > 1 )
    {
      cout.precision(6);

      cout << " Return value = new volume = "
           <<  (targetPhysical==0 ? G4String("None") :
		                    targetPhysical->GetName() )  << endl;
    }
#endif

#ifdef DEBUG_HIST
  G4cerr << "Upon exiting LocateGlobalPointAndSetup " << endl;
  G4cerr << "  History = " << endl << fHistory << endl << endl;
#endif
  return targetPhysical;	
}

// Compute the next geometric Step: Intersections with current
// mother and `daughter' volumes.
//
// NOTE:
//
// Flags on entry:
//
// fValidExitNormal  - Normal of exited volume is valid (convex, not a 
//                     coincident boundary)
// fExitNormal       - Surface normal of exited volume
// fExiting          - True if have exited solid
//
// fBlockedPhysicalVolume - Ptr to exited volume (or 0)
// fBlockedReplicaNo - Replication no of exited volume
// fLastStepWasZero  - True if last Step size was zero.
//
// Flags on exit:
// fValidExitNormal  - True if surface normal of exited volume is valid
// fExitNormal       - Surface normal of exited volume rotated to mothers
//                    reference system
// fExiting          - True if exiting mother
// fEntering         - True if entering `daughter' volume (or replica)
// fBlockedPhysicalVolume - Ptr to candidate (entered) volume
// fBlockedReplicaNo - Replication no of candidate (entered) volume
// fLastStepWasZero  - True if this Step size was zero.

G4double G4Navigator::ComputeStep(const G4ThreeVector &pGlobalpoint,
				  const G4ThreeVector &pDirection,
				  const G4double pCurrentProposedStepLength,
				  G4double	&pNewSafety)
{
  G4double Step;
  G4ThreeVector localDirection=ComputeLocalAxis(pDirection);
  G4VPhysicalVolume  *motherPhysical=fHistory.GetTopVolume();
  G4LogicalVolume *motherLogical=motherPhysical->GetLogicalVolume();

#ifdef G4VERBOSE
  cout.precision(8);
  if( fVerbose > 1 ) 
    {
      cout << "*** G4Navigator::ComputeStep: ***" << endl; 
      cout.precision(8);
      cout << " I was called with the following arguments: " << endl
	   << " Globalpoint = " << setw(25) << pGlobalpoint  << endl
	   << " Direction   = " << setw(25) << pDirection    << endl
	   << " ProposedStepLength= " << pCurrentProposedStepLength << endl;
       //       << " = " << << endl
    }

  if( fVerbose > 2 )
    {
      // cout.precision(3);
      cout << " Upon entering my state is: " << endl;
      PrintState();
    }
#endif

  G4ThreeVector newLocalPoint =ComputeLocalPoint(pGlobalpoint);
  if( newLocalPoint != fLastLocatedPointLocal )
    {
      // Check whether the relocation is within safety
      //
      G4ThreeVector oldLocalPoint= fLastLocatedPointLocal;
      G4double moveLenSq= (newLocalPoint-oldLocalPoint).mag2();

      if (moveLenSq >= kCarTolerance*kCarTolerance){
        //
	//  The following checks only make sense if the move is larger 
	//  than the tolerance.
	// 
	G4ThreeVector OriginalGlobalpoint;
	OriginalGlobalpoint = fHistory.GetTopTransform().Inverse()
	  .TransformPoint(fLastLocatedPointLocal); 

        G4double shiftOriginSafSq= (fPreviousSftOrigin-pGlobalpoint).mag2();

#if 0
        // Reset point before computing safety 
        LocateGlobalPointWithinVolume(OriginalGlobalpoint); 
	G4double safety= ComputeSafety(OriginalGlobalpoint);
	if( moveLenSq >= sqr(safety) ){
	   G4double moveLen=sqrt(moveLenSq);
	   if( moveLen > safety + kCarTolerance ){
	      G4cerr << " ERROR in G4Navigator::ComputeStep: " << endl
		   << "The Step's starting point has moved " << moveLen 
		   << " since the last call to one of the Locate methods " << endl
		   << " which is more than the current safety=" << safety << endl;
	   }else{
 	      G4cerr << " Warning in G4Navigator::ComputeStep: " << endl
		   << "The Step's starting point has moved " << moveLen 
		   << " which is equal to the current safety. " << endl;
	   }
	}
	G4double safetyPlus = safety + kCarTolerance;
	assert( moveLenSq <= sqr(safetyPlus) );
#endif 	
	if(  shiftOriginSafSq >= sqr(fPreviousSafety) ){
	   G4double shiftOrigin=sqrt(shiftOriginSafSq);
	   if( shiftOrigin > fPreviousSafety + kCarTolerance ){
	      G4cerr << " ERROR in G4Navigator::ComputeStep: " << endl
		   << "The Step's starting point has moved " << sqrt(moveLenSq)
		   << " since the last call to one of the Locate methods " << endl
		   << " This has resulted in moving " << shiftOrigin
		   << " from the last point at which the safety was calculated "
		   << endl 
		   << " which is more than the computed safety= " 
		   << fPreviousSafety << "at that point." << endl;
	   }
#ifdef DEBUG
	   else
	   {
 	      G4cerr << " Warning in G4Navigator::ComputeStep: " << endl
		   << "The Step's starting point has moved " << sqrt(moveLenSq)
		   << " which has taken it to the limit of the current safety. "
		   << endl;
	   }
#endif
	}
	G4double safetyPlus = fPreviousSafety+ kCarTolerance;
	assert( shiftOriginSafSq <= sqr(safetyPlus) );

	// Relocate the point within the same volume
	//
	LocateGlobalPointWithinVolume( pGlobalpoint );
      }
    }

  if (fHistory.GetTopVolumeType()!=kReplica)
    {
      switch(CharacteriseDaughters(motherLogical))
	{
	case kNormal:
	  if (motherLogical->GetVoxelHeader())
	    {
	      Step=fvoxelNav.ComputeStep(fLastLocatedPointLocal,
					 localDirection,
					 pCurrentProposedStepLength,
					 pNewSafety,
					 fHistory,
					 fValidExitNormal,
					 fExitNormal,
					 fExiting,
					 fEntering,
					 &fBlockedPhysicalVolume,
					 fBlockedReplicaNo);
	      
	    }
	  else
	    {
	      Step=fnormalNav.ComputeStep(fLastLocatedPointLocal,
					  localDirection,
					  pCurrentProposedStepLength,
					  pNewSafety,
					  fHistory,
					  fValidExitNormal,
					  fExitNormal,
					  fExiting,
					  fEntering,
					  &fBlockedPhysicalVolume,
					  fBlockedReplicaNo);
	    }
	  break;
	case kParameterised:
	  Step=fparamNav.ComputeStep(fLastLocatedPointLocal,
				     localDirection,
				     pCurrentProposedStepLength,
				     pNewSafety,
				     fHistory,
				     fValidExitNormal,
				     fExitNormal,
				     fExiting,
				     fEntering,
				     &fBlockedPhysicalVolume,
				     fBlockedReplicaNo);
	  break;
	case kReplica:
	  G4Exception("Logic Error in G4Navigator::ComputeStep()");
	  break;
	}
    }
  else
    {
      // In the case of a replica, 
      //   it must handles the exiting edge/corner problem by itself
      G4bool exitingReplica= fExitedMother;
      Step=freplicaNav.ComputeStep(pGlobalpoint,
				   pDirection,
				   fLastLocatedPointLocal,
				   localDirection,
				   pCurrentProposedStepLength,
				   pNewSafety,
				   fHistory,
				   fValidExitNormal,
				   fExitNormal,
				   exitingReplica,
				   fEntering,
				   &fBlockedPhysicalVolume,
				   fBlockedReplicaNo);
      // still ok to set it ??
      fExiting= exitingReplica;
     }

  if( (Step == pCurrentProposedStepLength) && (!fExiting) && (!fEntering) )
    {
      // This is Step is not really limited by the geometry.
      // The Navigator is obliged to return "infinity"
      Step = kInfinity;
    }

  // Remember last safety origin & value.
  fPreviousSftOrigin= pGlobalpoint;
  fPreviousSafety= pNewSafety; 

  fLocatedOnEdge= fLastStepWasZero && (Step==0);  // Edge if two consecutive
                                                  // steps are zero, because
                   //  at least two candidate volumes must have been checked

  fLastStepWasZero= (Step==0);
  fEnteredDaughter=fEntering;   // I expect to enter a volume in this Step
  fExitedMother=fExiting;

  if(fExiting && !fValidExitNormal)
     {
       // We must calculate the normal anyway (in order to have it if requested)
       G4ThreeVector FinalPoint=  fLastLocatedPointLocal + localDirection*Step;
       fExitNormal= motherLogical->GetSolid()->SurfaceNormal(FinalPoint);
     }

#ifdef G4VERBOSE
  if( fVerbose > 1 ) 
    {
      cout << " Upon exiting my state is: " << endl;
      PrintState();
    }
#endif

  return Step;
}

G4VPhysicalVolume* G4Navigator::LocateGlobalPointAndSetup(const G4ThreeVector &p,
							  const G4TouchableHistory &h)
{
  fHistory=*h.GetHistory();
  SetupHierarchy();
  return LocateGlobalPointAndSetup(p, 0);
}

G4ThreeVector G4Navigator::NetTranslation() const
{
  G4AffineTransform tf(fHistory.GetTopTransform().Inverse());
  return tf.NetTranslation();
}

G4RotationMatrix G4Navigator::NetRotation() const
{
  G4AffineTransform tf(fHistory.GetTopTransform().Inverse());
  return tf.NetRotation();
}

G4GRSVolume* G4Navigator::CreateGRSVolume() const
{
  G4AffineTransform tf(fHistory.GetTopTransform().Inverse());
  return new G4GRSVolume(fHistory.GetTopVolume(),
			 tf.NetRotation(),
			 tf.NetTranslation());
}

G4GRSSolid* G4Navigator::CreateGRSSolid() const
{
  G4AffineTransform tf(fHistory.GetTopTransform().Inverse());
  return new G4GRSSolid(fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid(),
			tf.NetRotation(),
			tf.NetTranslation());
			
}

G4TouchableHistory* G4Navigator::CreateTouchableHistory() const
{
  return new G4TouchableHistory(fHistory);
}

// Renavigate & reset hierarchy described by current history
// o Reset volumes
// o Recompute transforms and/or solids of replicated/parameterised vols
void G4Navigator::SetupHierarchy()
{
  G4int i;
  const G4int cdepth=fHistory.GetDepth();
  G4VPhysicalVolume *mother,*current;
  G4VSolid *pSolid;
  G4VPVParameterisation *pParam;

  mother=fHistory.GetVolume(0);
  for (i=1;i<=cdepth;i++)
  {
    current=fHistory.GetVolume(i);
    switch (fHistory.GetVolumeType(i))
      {
      case kNormal:
	break;
      case kReplica:
	freplicaNav.ComputeTransformation(fHistory.GetReplicaNo(i),
					  current);
	break;
      case kParameterised:
        G4int replicaNo;
	// pSolid=current->GetLogicalVolume()->GetSolid();
	pParam=current->GetParameterisation();
        replicaNo= fHistory.GetReplicaNo(i);
	pSolid= pParam->ComputeSolid(replicaNo, current);
	// Set up dimensions & transform in solid/physical volume
	pSolid->ComputeDimensions(pParam, replicaNo, current);
	pParam->ComputeTransformation(replicaNo, current);
	
	// Set up the correct solid and material in Logical Volume
	G4LogicalVolume *pLogical;
	pLogical= current->GetLogicalVolume();
	pLogical->SetSolid( pSolid );
	pLogical->SetMaterial( pParam->ComputeMaterial(replicaNo, 
					               current));
	break;
      }
    current->Setup(mother);
    mother=current;
  }
}

ostream& operator << (ostream &os,const G4Navigator &n)
{

  os << "Current History: " << endl << n.fHistory;
  return os;
}

// Return global to local transformation 
const G4AffineTransform G4Navigator::GetLocalToGlobalTransform() const
{
  G4AffineTransform  tempTransform;
  tempTransform= fHistory.GetTopTransform().Inverse(); 
  return  tempTransform;
}

//  Obtain the Normal vector to a surface (in local coordinates)
//   pointing out of previous volume and into current volume
//    
G4ThreeVector  G4Navigator::GetLocalExitNormal(G4bool* valid)
{
  G4ThreeVector ExitNormal(0.,0.,0.);

  if( fExitedMother ){
     ExitNormal=fExitNormal;
     *valid = true;
     
  }else if (EnteredDaughterVolume()) {
     ExitNormal= -(fHistory.GetTopVolume()->GetLogicalVolume()
			 ->GetSolid()->SurfaceNormal(fLastLocatedPointLocal));
     *valid = true;
  }else{
     // We are not at a boundary.
     // ExitNormal remains (0,0,0) 
     *valid = false;
  }

  return ExitNormal;
}

//   It assumes that it assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.

G4double G4Navigator::ComputeSafety(const G4ThreeVector &pGlobalpoint,
				  const G4double pMaxLength)
                                               // A sort of MaximumLength ... ?
{
  G4double newSafety=0.0;

#ifdef G4VERBOSE
  if( fVerbose > 0 ) 
    {
      cout << "*** G4Navigator::ComputeSafety: ***" << endl; 
      cout.precision(8);
      cout << " I was called with the following arguments: " << endl
	   << " Globalpoint = " << pGlobalpoint << endl;
      //       cout << "  pMaxLength  = " << pMaxLength  << endl;

      cout << " Upon entering my state is: " << endl;
      PrintState();
    }
#endif

  // Pseudo-relocate to this point (updates voxel information only).
  LocateGlobalPointWithinVolume( pGlobalpoint );

  if( ! (fEnteredDaughter || fExitedMother ) )
    {
        G4VPhysicalVolume  *motherPhysical=fHistory.GetTopVolume();
        G4LogicalVolume *motherLogical=motherPhysical->GetLogicalVolume();

        G4ThreeVector localPoint= ComputeLocalPoint(pGlobalpoint);
	if (fHistory.GetTopVolumeType()!=kReplica)
	  {
	    switch(CharacteriseDaughters(motherLogical))
	      {
	      case kNormal:
		if (motherLogical->GetVoxelHeader())
		  {
		    newSafety=fvoxelNav.ComputeSafety(localPoint,
					       fHistory,	      
					       pMaxLength);
		  }
		else
		  {
	    
		    newSafety=fnormalNav.ComputeSafety(localPoint,
						       fHistory,	      
						       pMaxLength);

		  }
		break;
	      case kParameterised:

		newSafety=fparamNav.ComputeSafety(localPoint,
					          fHistory,	      
					          pMaxLength);
		break;
	      case kReplica:
		G4Exception("Logic Error in G4Navigator::ComputeSafety()");
		break;
	      }
	  }
	else
	  {
	    newSafety=freplicaNav.ComputeSafety(pGlobalpoint,
						localPoint,
					        fHistory,	      
					        pMaxLength);                
	  }
    }

  // Remember last safety origin & value.
  fPreviousSftOrigin= pGlobalpoint;
  fPreviousSafety= newSafety; 

#ifdef G4VERBOSE
  if( fVerbose > 1 ) 
    {
      cout.precision(8);
      cout << " Upon exiting my state is: " << endl;
      PrintState();
      cout << "  and I return a value of Safety = " << newSafety << endl;
    }
#endif


  return newSafety;
}

G4bool  G4Navigator::EnteredDaughterVolume()
{
  return fEnteredDaughter;
}

// G4bool        G4Navigator::ExitedVolume()
// {
//   return fExitedCurrent;
// }


void  G4Navigator::PrintState()
{
  if( fVerbose >= 4 )
    {
      cout.precision(3);
      cout << " Upon exiting my state is: " << endl;
      cout << "  ValidExitNormal= " << fValidExitNormal << endl
	   << "  ExitNormal     = " << fExitNormal      << endl
	   << "  Exiting        = " << fExiting         << endl
	   << "  Entering       = " << fEntering        << endl
	   << "  BlockedPhysicalVolume= " <<  (fBlockedPhysicalVolume==0 ? G4String("None") :
					       fBlockedPhysicalVolume->GetName() )              << endl
	   << "  BlockedReplicaNo     = " <<  fBlockedReplicaNo       << endl
	   << "  LastStepWasZero      = " <<   fLastStepWasZero       << endl
	   << endl;   
    }
  if( ( 1 < fVerbose) && (fVerbose < 4) )
    {
      cout.precision(3);
      cout << setw(18) << " ExitNormal "  << " "     
	   << setw( 5) << " Valid "       << " "     
	   << setw( 9) << " Exiting "     << " "      
	   << setw( 9) << " Entering"     << " " 
	   << setw(15) << " Blocked:Volume "  << " "   
	   << setw( 9) << " ReplicaNo"        << " "  
	   << setw( 8) << " LastStepZero  "   << " "   
	   << endl;   
      cout << setw(18)  << fExitNormal       << " "
	   << setw( 5)  << fValidExitNormal  << " "   
	   << setw( 9)  << fExiting          << " "
	   << setw( 9)  << fEntering         << " "
	   << setw(15)  << (fBlockedPhysicalVolume==0 ? G4String("None") :
				    fBlockedPhysicalVolume->GetName() )   << " "   
	   << setw( 9)  << fBlockedReplicaNo  << " "
	   << setw( 8)  << fLastStepWasZero   << " "
	   << endl;   
    }
  if( fVerbose > 2 ) 
    {
      cout.precision(8);
      cout << " Current Localpoint = " << fLastLocatedPointLocal << endl;
      cout << " PreviousSftOrigin  = " << fPreviousSftOrigin << endl;
      cout << " PreviousSafety     = " <<  fPreviousSafety << endl; 
    }
}

void G4Navigator::LocateGlobalPointWithinVolume(const  G4ThreeVector& pGlobalpoint)
{
  // The new implementation of LocateGlobalPointWithinVolume 
  // 
  //   -> the state information of this Navigator and its subNavigators
  //       is updated in order to start the next step at pGlobalpoint
  //   -> no check is performed whether pGlobalpoint is inside the 
  //       original volume (this must be the case).
  //
  //  Note: a direction could be added to the arguments, to aid in
  //        future optional checking (via the Old code below). 
  //        [ This would be done only in verbose mode ]
  
   fLastLocatedPointLocal =ComputeLocalPoint(pGlobalpoint);

   // For the case of Voxel (or Parameterised) volume the respective 
   //   Navigator must be messaged to update its voxel information etc.o

   // Update the state of the Sub Navigators 
   //   - in particular any voxel information they store/cache
   //.
   G4VPhysicalVolume*  motherPhysical=fHistory.GetTopVolume();
   G4LogicalVolume*    motherLogical= motherPhysical->GetLogicalVolume();
   G4SmartVoxelHeader* pVoxelHeader=  motherLogical->GetVoxelHeader();

   G4ThreeVector localPoint= ComputeLocalPoint(pGlobalpoint);
   if (fHistory.GetTopVolumeType()!=kReplica)
     {
       switch(CharacteriseDaughters(motherLogical))
	 {
	 case kNormal:
	   if (pVoxelHeader)
	     {
	       fvoxelNav.VoxelLocate( pVoxelHeader, localPoint );
	     }
	     //  else { fnormalNav. nothing !? }
	   break;

	 case kParameterised:
	   // Resets state & returns voxel node 
	   fparamNav.VoxelLocate( pVoxelHeader, localPoint );
	   break;

	 case kReplica:
	   G4Exception("Logic Error in G4Navigator::LocateGlobalPointWithinVolume()");
	   break;
	 }
     }

#if 0
   else
     {
       // There is no state stored in G4ReplicaNavigation
       // freplicaNav.VoxelLocate( pVoxelHeader, localPoint );                
     }
#endif


#ifdef OLD_LOCATE
   //  An alternative implementation using LocateGlobalPointAndSetup.
   //  It can also be used to check the method's assumptions. 
   // 
   G4VPhysicalVolume  *pOldVol, *pNewVol;

   pOldVol= fHistory.GetTopVolume(); 
   pNewVol= LocateGlobalPointAndSetup(pGlobalpoint, 0); 
                                               // , G4ThreeVector(1.,0.,0.));

   if( pOldVol != pNewVol ){
      // This is abnormal behaviour.
      cerr << " ERROR in G4Navigator::LocateGlobalPointWithinVolume " << endl;
      cerr << "       A volume change has occured - this is not expected & illegal" << endl; 
      cerr << "         Old volume name = " << pOldVol->GetName() << endl;
      cerr << "         New volume name = " << pNewVol->GetName() << endl;

      G4VPhysicalVolume  *pNewVol2;
      pNewVol2= LocateGlobalPointAndSetup(pGlobalpoint, 0); 
                                                   //, G4ThreeVector(1.,0.,0.));
      cerr << "         Tried again & found volume= " << pNewVol2->GetName() << endl;

   }

   // Check that the new volume located is same as the old one.
   assert( pOldVol == pNewVol );
#endif

}

G4int G4Navigator::GetVerboseLevel()
{
   return fVerbose;
}

void  G4Navigator::SetVerboseLevel(G4int level)
{
   fVerbose=level;
}

