// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Navigator.cc,v 1.10 1999-12-15 17:51:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Navigator Implementation  Paul Kent July 95/96

#include "G4Navigator.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

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
  G4cerr << "Upon entering LocateGlobalPointAndSetup " << G4endl;
  G4cerr << "  History = " << G4endl << fHistory << G4endl << G4endl;
#endif
#ifdef G4VERBOSE
  if( fVerbose > 0 ) 
    {
      G4cout << "*** G4Navigator::LocateGlobalPointAndSetup: ***" << G4endl; 
      G4cout.precision(8);
      G4cout << " I was called with the following arguments: " << G4endl
	     << " Globalpoint = " << globalPoint << G4endl
	     << " relativeSearch   = " <<  relativeSearch  << G4endl;
       //       << " = " << << G4endl
      G4cout << " Upon entering my state is: " << G4endl;
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
	      //  - This stops it from exiting further volumes and cycling
	      //  - However ReplicaNavigator treats this case itself
	      if( fLocatedOnEdge 
                  && (VolumeType(fBlockedPhysicalVolume) != kReplica ))
                                                    // ( fLastStepWasZero )
		{ 
		  fExiting= false;
		}
	    }
	  else if (fEntering)
	    {
	      G4VPhysicalVolume *curPhysical=fHistory.GetTopVolume();
	      switch (VolumeType(fCandidatePhysicalVolume))
		{
		case kNormal:
		  fCandidatePhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fCandidatePhysicalVolume);
		  break;
		case kReplica:
		  freplicaNav.ComputeTransformation(fCandidateReplicaNo,
						    fCandidatePhysicalVolume);
		  fCandidatePhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fCandidatePhysicalVolume,
				    kReplica,
				    fCandidateReplicaNo);
		  fCandidatePhysicalVolume->SetCopyNo(fCandidateReplicaNo);
		  
		  break;
		case kParameterised:
		  G4VSolid *pSolid;
		  // G4VSolid *pSolid=fCandidatePhysicalVolume->
		  // GetLogicalVolume()-> GetSolid();
		  G4VPVParameterisation *pParam=fCandidatePhysicalVolume->
		    GetParameterisation();
		  pSolid= pParam->ComputeSolid(fCandidateReplicaNo,
						fCandidatePhysicalVolume);
		  pSolid->ComputeDimensions(pParam,
					    fCandidateReplicaNo,
					    fCandidatePhysicalVolume);
		  pParam->ComputeTransformation(fCandidateReplicaNo,
						fCandidatePhysicalVolume);
		  fCandidatePhysicalVolume->Setup(curPhysical);
		  fHistory.NewLevel(fCandidatePhysicalVolume,
				    kParameterised,
				    fCandidateReplicaNo);
		  fCandidatePhysicalVolume->SetCopyNo(fCandidateReplicaNo);
		  // Set the correct solid and material in Logical Volume
		  G4LogicalVolume *pLogical;
		  pLogical= fCandidatePhysicalVolume->GetLogicalVolume();
		  pLogical->SetSolid( pSolid );
		  pLogical->SetMaterial( 
			     pParam->ComputeMaterial(fCandidateReplicaNo, 
						     fCandidatePhysicalVolume));
		  break;
		}
	      fEntering=false;
	      fCandidatePhysicalVolume=0;
	      localPoint=fHistory.GetTopTransform().TransformPoint(globalPoint);
	      notKnownContained=false;
	    }
	}
      else
	{
	  fBlockedPhysicalVolume=0;
	  fCandidatePhysicalVolume=0;
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
      G4cout.precision(6);

      G4String curPhysVol_Name("None");
      if (targetPhysical!=0)
	 curPhysVol_Name= targetPhysical->GetName();
      G4cout << " Return value = new volume = "
	     << curPhysVol_Name  << G4endl;
    }
#endif

#ifdef DEBUG_HIST
  G4cerr << "Upon exiting LocateGlobalPointAndSetup " << G4endl;
  G4cerr << "  History = " << G4endl << fHistory << G4endl << G4endl;
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
// fCandidatePhysicalVolume - Ptr to candidate (entered) volume
// fCandidateReplicaNo - Replication no of candidate (entered) volume
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
  G4cout.precision(8);
  if( fVerbose > 1 ) 
    {
      G4cout << "*** G4Navigator::ComputeStep: ***" << G4endl; 
      G4cout.precision(8);
      G4cout << " I was called with the following arguments: " << G4endl
	   << " Globalpoint = " << G4std::setw(25) << pGlobalpoint  << G4endl
	   << " Direction   = " << G4std::setw(25) << pDirection    << G4endl
	   << " ProposedStepLength= " << pCurrentProposedStepLength << G4endl;
       //       << " = " << << G4endl
    }

  if( fVerbose > 2 )
    {
      // G4cout.precision(3);
      G4cout << " Upon entering my state is: " << G4endl;
      PrintState();
    }
#endif

  static G4double fAccuracyForWarning=   kCarTolerance,
                  fAccuracyForException= 1000*kCarTolerance;

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

	//  Check that the starting point of this step is 
	//   within the isotropic safety sphere of the last point
	//   to a accuracy/precision  given by 
	//       fAccuracyForWarning
	//   If so give warning.  If it fails by more than
	//       fAccuracyForException
        //   exit with error.
  
	if(  shiftOriginSafSq >= sqr(fPreviousSafety) ){
	   G4double shiftOrigin=sqrt(shiftOriginSafSq);
	   G4double diffShiftSaf=  shiftOrigin - fPreviousSafety;
	   G4bool isError; 

	   if( diffShiftSaf > fAccuracyForWarning ){
              isError = ( diffShiftSaf >= fAccuracyForException );
	      G4cerr.precision(10);
              if ( isError )
		G4cerr << "Accuracy ERROR found in G4Navigator::ComputeStep: " << G4endl;
              else
		G4cerr << "Warning G4Navigator::ComputeStep found slightly inaccurate position:" << G4endl;

	      G4cerr << "     The Step's starting point has moved " 
		     << sqrt(moveLenSq)/mm << " mm " << G4endl
		     << "     since the last call to a Locate method." << G4endl;
	      G4cerr << "     This has resulted in moving " 
		     << shiftOrigin/mm << " mm " 
		     << " from the last point at which the safety " 
		     << "     was calculated " << G4endl;
	      G4cerr << "     which is more than the computed safety= " 
		     << fPreviousSafety/mm << " mm  at that point." << G4endl;
	      G4cerr << "     This difference is " 
		     << diffShiftSaf /mm << " mm." << G4endl;

#ifdef G4VERBOSE
              static G4int warnNow= 0;
	      if( ((++warnNow % 100) == 1) ) {   //  || (warnNow < 4) ){ 
	          G4cerr << "  This problem can be due to either " << G4endl;
		  G4cerr << "    - a process that has proposed a displacement"
			 <<     " larger than the current safety , or" << G4endl;
		  G4cerr << "    - inaccuracy in the computation of the safety" << G4endl;
		  G4cerr << "    - if you are using a magnetic field, a known conflict about the safety exists in this case."
			 << G4endl;
		  G4cerr << "  We suggest that you " << G4endl
                     <<  "   - find i) what particle is being tracked, and "
			 << " ii) through what part of your geometry " << G4endl
		     <<  "      for example by reruning this event with " << G4endl
		     <<  "         /tracking/verbose 1 "  << G4endl
                     << "    - check which processes you declare for this particle"
			 << " (and look at non-standard ones) "  << G4endl
                     <<  "   - if possible create a detailed logfile "
		         << " of this event using:" << G4endl
                     <<  "         /tracking/verbose 6 "
                     << G4endl;
	      }
              // G4cerr << "    -  ." << G4endl;
#endif 
	   }
#ifdef DEBUG
	   else
	   {
 	      G4cerr << " Warning in G4Navigator::ComputeStep: " << G4endl
		   << "The Step's starting point has moved " << sqrt(moveLenSq)
		   << " which has taken it to the limit of the current safety. "
		   << G4endl;
	   }
#endif
	}
	G4double safetyPlus = fPreviousSafety+ fAccuracyForException;
	assert( shiftOriginSafSq <= sqr(safetyPlus) );

	// Relocate the point within the same volume
	//
	LocateGlobalPointWithinVolume( pGlobalpoint );
      }
    }

  G4VPhysicalVolume *lBlockedOrCandidatePhysicalVolume;
  G4int lBlockedOrCandidateReplicaNo;

  // We need to give the "blocked volume" as input to the subNavigators
  //   and receive the "Candidate volume"
  lBlockedOrCandidatePhysicalVolume = fBlockedPhysicalVolume;
  lBlockedOrCandidateReplicaNo      = fBlockedReplicaNo;

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
					 &lBlockedOrCandidatePhysicalVolume,
					 lBlockedOrCandidateReplicaNo);
	      
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
					  &lBlockedOrCandidatePhysicalVolume,
					  lBlockedOrCandidateReplicaNo);
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
				     &lBlockedOrCandidatePhysicalVolume,
				     lBlockedOrCandidateReplicaNo);
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
				   &lBlockedOrCandidatePhysicalVolume,
				   lBlockedOrCandidateReplicaNo);
      // still ok to set it ??
      fExiting= exitingReplica;
     }

  // We receive the "Candidate volume"
  fCandidatePhysicalVolume= lBlockedOrCandidatePhysicalVolume;
  fCandidateReplicaNo     = lBlockedOrCandidateReplicaNo;  

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
      G4cout << " Upon exiting my state is: " << G4endl;
      PrintState();
      G4cout << "--> am returning a Step= " << G4std::setw(6) << Step << G4endl 
	     << G4endl;
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

G4std::ostream& operator << (G4std::ostream &os,const G4Navigator &n)
{

  os << "Current History: " << G4endl << n.fHistory;
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
      G4cout << "*** G4Navigator::ComputeSafety: ***" << G4endl; 
      G4cout.precision(5);
      G4cout << " I was called with the following arguments: " << G4endl
	     << " Globalpoint = " 
	     << G4std::setw(24) << pGlobalpoint << G4endl;
      //       G4cout << "  pMaxLength  = " << pMaxLength  << G4endl;

      G4cout << " Upon entering my state is: " << G4endl;
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
      G4cout.precision(8);
      G4cout << " Upon exiting my state is: " << G4endl;
      PrintState();
      G4cout << "  and I return a value of Safety = " << newSafety << G4endl;
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
      G4cout.precision(3);
      G4cout << " Upon exiting my state is: " << G4endl;
      G4cout << "  ValidExitNormal= " << fValidExitNormal << G4endl
	   << "  ExitNormal     = " << fExitNormal      << G4endl
	   << "  Exiting        = " << fExiting         << G4endl
	   << "  Entering       = " << fEntering        << G4endl
	   << "  BlockedPhysicalVolume= " ;
      if (fBlockedPhysicalVolume==0 )
	 G4cout << "None";
      else
 	 G4cout << fBlockedPhysicalVolume->GetName();
      G4cout << G4endl
	   << "  BlockedReplicaNo     = " <<  fBlockedReplicaNo       << G4endl
	   << "  LastStepWasZero      = " <<   fLastStepWasZero       << G4endl
	   << G4endl;   
    }
  if( ( 1 < fVerbose) && (fVerbose < 4) )
    {
      G4cout.precision(3);
      G4cout << G4std::setw(15) << "   ExitNormal  "  << " "     
	   << G4std::setw( 5) << "Valid"       << " "     
 	   << G4std::setw(10) << "Exit/Enter"     << " "        //  Values use 7
	   << G4std::setw(15) << "Blocked:Volume "  << " "   
	   << G4std::setw( 5) << "Repl#"        << " "  
	   << G4std::setw( 7) << "lStep=0"   << " "   
	   << G4std::setw(12) << "Candidate:Volume "  << " "   
	   << G4endl;   

      G4cout << G4std::setw(10)  << fExitNormal       << " "
	   << G4std::setw( 4)  << fValidExitNormal  << " "   
	   << G4std::setw( 3)  << fExiting          << " "
	   << G4std::setw( 2)  << fEntering         << " ";
      if (fBlockedPhysicalVolume==0 )
	 G4cout << G4std::setw(12) << "None";
      else
 	 G4cout << G4std::setw(12)<< fBlockedPhysicalVolume->GetName();
      G4cout << G4std::setw( 8) << fBlockedPhysicalVolume << " "; // Hex would be better
      G4cout << G4std::setw( 6)  << fBlockedReplicaNo  << " "
	     << G4std::setw( 4)  << fLastStepWasZero   << " "; 
      if (fCandidatePhysicalVolume==0 )
	 G4cout << G4std::setw(12) << "None";
      else
 	 G4cout << G4std::setw(12)<< fCandidatePhysicalVolume->GetName();
      G4cout << G4std::setw( 4)  << fCandidateReplicaNo  << " ";
      G4cout << G4std::setw(12)  << fCandidatePhysicalVolume << " "
	     << G4endl;   
    }
  if( fVerbose > 2 ) 
    {
      G4cout.precision(8);
      G4cout << " Current Localpoint = " << fLastLocatedPointLocal << G4endl;
      G4cout << " PreviousSftOrigin  = " << fPreviousSftOrigin << G4endl;
      G4cout << " PreviousSafety     = " <<  fPreviousSafety << G4endl; 
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

   else
     {
       // There is no state stored in G4ReplicaNavigation
       // freplicaNav.VoxelLocate( pVoxelHeader, localPoint );                
     }

   // Since we are doing a (fast) relocation,  invalidate 
   //  any old values of fEntering, fExiting - which could no longer hold!
   //                            July 20, 1999
   fEntering = false;
   fExiting  = false;

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
      G4cerr << " ERROR in G4Navigator::LocateGlobalPointWithinVolume " << G4endl;
      G4cerr << "       A volume change has occured - this is not expected & illegal" << G4endl; 
      G4cerr << "         Old volume name = " << pOldVol->GetName() << G4endl;
      G4cerr << "         New volume name = " << pNewVol->GetName() << G4endl;

      G4VPhysicalVolume  *pNewVol2;
      pNewVol2= LocateGlobalPointAndSetup(pGlobalpoint, 0); 
                                                   //, G4ThreeVector(1.,0.,0.));
      G4cerr << "         Tried again & found volume= " << pNewVol2->GetName() << G4endl;

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

