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
// class G4ParameterisedNavigation Implementation
//
// Initial Author: P.Kent, 1996
// Revisions:
//  J. Apostolakis 24 Nov 2005, Revised/fixed treatment of nested params
//  J. Apostolakis  4 Feb 2005, Reintroducting multi-level parameterisation
//                              for materials only - see note 1 below
//  G. Cosmo       11 Mar 2004, Added Check mode 
//  G. Cosmo       15 May 2002, Extended to 3-d voxelisation, made subclass
//  J. Apostolakis  5 Mar 1998, Enabled parameterisation of mat & solid type
// --------------------------------------------------------------------

// Note 1: Design/implementation note for extensions - JAp, March 1st, 2005
// We cannot make the solid, dimensions and transformation dependent on
// parent because the voxelisation will not have access to this. 
// So the following can NOT be done:
//   sampleSolid = curParam->ComputeSolid(num, curPhysical, pParentTouch);
//   sampleSolid->ComputeDimensions(curParam, num, curPhysical, pParentTouch);
//   curParam->ComputeTransformation(num, curPhysical, pParentTouch);

#include "G4ParameterisedNavigation.hh"
#include "G4TouchableHistory.hh"
#include "G4VNestedParameterisation.hh"

#include "G4AuxiliaryNavServices.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4ParameterisedNavigation::G4ParameterisedNavigation()
{
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4ParameterisedNavigation::~G4ParameterisedNavigation()
{
}

// ***************************************************************************
// ComputeStep
// ***************************************************************************
//
G4double G4ParameterisedNavigation::
                    ComputeStep(const G4ThreeVector& localPoint,
                                const G4ThreeVector& localDirection,
                                const G4double currentProposedStepLength,
                                      G4double& newSafety,
                                      G4NavigationHistory& history,
                                      G4bool& validExitNormal,
                                      G4ThreeVector& exitNormal,
                                      G4bool& exiting,
                                      G4bool& entering,
                                      G4VPhysicalVolume *(*pBlockedPhysical),
                                      G4int& blockedReplicaNo)
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical;
  G4VPVParameterisation *sampleParam;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid, *sampleSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, ourSafety;
  G4double motherSafety, motherStep = DBL_MAX;
  G4bool motherValidExitNormal = false;
  G4ThreeVector motherExitNormal;
  
  G4int sampleNo;

  G4bool initialNode, noStep;
  G4SmartVoxelNode *curVoxelNode;
  G4long curNoVolumes, contentNo;
  G4double voxelSafety;

  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;

  motherPhysical = history.GetTopVolume();
  motherLogical = motherPhysical->GetLogicalVolume();
  motherSolid = motherLogical->GetSolid();

  //
  // Compute mother safety
  //

  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety;              // Working isotropic safety

#ifdef G4VERBOSE
  if ( fCheck )
  {
    if( motherSafety < 0.0 )
    {
      motherSolid->DumpInfo();
      std::ostringstream message;
      message << "Negative Safety In Voxel Navigation !" << G4endl
              << "        Current solid " << motherSolid->GetName()
              << " gave negative safety: " << motherSafety << G4endl
              << "        for the current (local) point " << localPoint;
      G4Exception("G4ParameterisedNavigation::ComputeStep()",
                  "GeomNav0003", FatalException, message); 
    }
    if( motherSolid->Inside(localPoint) == kOutside )
    { 
      std::ostringstream message;
      message << "Point is outside Current Volume !" << G4endl
              << "          Point " << localPoint
              << " is outside current volume " << motherPhysical->GetName()
              << G4endl;
      G4double estDistToSolid = motherSolid->DistanceToIn(localPoint); 
      G4cout << "          Estimated isotropic distance to solid (distToIn)= " 
             << estDistToSolid;
      if( estDistToSolid > 100.0 * motherSolid->GetTolerance() )
      {
        motherSolid->DumpInfo();
        G4Exception("G4ParameterisedNavigation::ComputeStep()",
                    "GeomNav0003", FatalException, message,
                    "Point is far outside Current Volume !"); 
      }
      else
      {
        G4Exception("G4ParameterisedNavigation::ComputeStep()",
                    "GeomNav1002", JustWarning, message,
                    "Point is a little outside Current Volume.");
      }
    }

    // Compute early:
    //  a) to check whether point is (wrongly) outside
    //               (signaled if step < 0 or step == kInfinity )
    //  b) to check value against answer of daughters!
    //
    motherStep = motherSolid->DistanceToOut(localPoint,
                                            localDirection,
                                            true,
                                           &motherValidExitNormal,
                                           &motherExitNormal);
  
    if( (motherStep >= kInfinity) || (motherStep < 0.0) )
    {
      // Error - indication of being outside solid !!
      //
      fLogger->ReportOutsideMother(localPoint, localDirection, motherPhysical);

      ourStep = motherStep = 0.0;
      exiting = true;
      entering = false;
    
      // If we are outside the solid does the normal make sense?
      validExitNormal = motherValidExitNormal;
      exitNormal = motherExitNormal;
    
      *pBlockedPhysical = nullptr; // or motherPhysical ?
      blockedReplicaNo = 0;  // or motherReplicaNumber ?
    
      newSafety = 0.0;
      return ourStep;
    }
  }
#endif

  initialNode = true;
  noStep = true;

  // By definition, the parameterised volume is the first
  // (and only) daughter of the mother volume
  //
  samplePhysical = motherLogical->GetDaughter(0);
  samplePhysical->GetReplicationData(axis,nReplicas,width,offset,consuming);
  fBList.Enlarge(nReplicas);
  fBList.Reset();

  // Exiting normal optimisation
  //
  if (exiting && (*pBlockedPhysical==samplePhysical) && validExitNormal)
  {
    if (localDirection.dot(exitNormal)>=kMinExitingNormalCosine)
    {
      // Block exited daughter replica; Must be on boundary => zero safety
      //
      fBList.BlockVolume(blockedReplicaNo);
      ourSafety = 0;
    }
  }
  exiting = false;
  entering = false;

  sampleParam = samplePhysical->GetParameterisation();

  // Loop over voxels & compute daughter safeties & intersections

  do
  {
    curVoxelNode = fVoxelNode;
    curNoVolumes = curVoxelNode->GetNoContained();

    for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
    {
      sampleNo = curVoxelNode->GetVolume((G4int)contentNo);
      if ( !fBList.IsBlocked(sampleNo) )
      {
        fBList.BlockVolume(sampleNo);

        // Call virtual methods, and copy information if needed
        //
        sampleSolid = IdentifyAndPlaceSolid( sampleNo, samplePhysical,
                                             sampleParam ); 

        G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                   samplePhysical->GetTranslation());
        sampleTf.Invert();
        const G4ThreeVector samplePoint = sampleTf.TransformPoint(localPoint);
        const G4double sampleSafety = sampleSolid->DistanceToIn(samplePoint);
        if ( sampleSafety<ourSafety )
        {
          ourSafety = sampleSafety;
        }
        if ( sampleSafety<=ourStep )
        {
          sampleDirection = sampleTf.TransformAxis(localDirection);
          G4double sampleStep =
                   sampleSolid->DistanceToIn(samplePoint, sampleDirection);
          if ( sampleStep<=ourStep )
          {
            ourStep = sampleStep;
            entering = true;
            exiting = false;
            *pBlockedPhysical = samplePhysical;
            blockedReplicaNo = sampleNo;
#ifdef G4VERBOSE
              // Check to see that the resulting point is indeed in/on volume.
              // This check could eventually be made only for successful
              // candidate.

              if ( ( fCheck ) && ( sampleStep < kInfinity ) )
              {
                G4ThreeVector intersectionPoint;
                intersectionPoint = samplePoint + sampleStep * sampleDirection;
                EInside insideIntPt = sampleSolid->Inside(intersectionPoint); 
                if( insideIntPt != kSurface )
                {
                  G4long oldcoutPrec = G4cout.precision(16); 
                  std::ostringstream message;
                  message << "Navigator gets conflicting response from Solid."
                          << G4endl
                          << "          Inaccurate solid DistanceToIn"
                          << " for solid " << sampleSolid->GetName() << G4endl
                          << "          Solid gave DistanceToIn = "
                          << sampleStep << " yet returns " ;
                  if( insideIntPt == kInside )
                    message << "-kInside-"; 
                  else if( insideIntPt == kOutside )
                    message << "-kOutside-";
                  else
                    message << "-kSurface-"; 
                  message << " for this point !" << G4endl
                          << "          Point = " << intersectionPoint
                          << G4endl;
                  if ( insideIntPt != kInside )
                    message << "        DistanceToIn(p) = " 
                            << sampleSolid->DistanceToIn(intersectionPoint);
                  if ( insideIntPt != kOutside ) 
                    message << "        DistanceToOut(p) = " 
                            << sampleSolid->DistanceToOut(intersectionPoint);
                  G4Exception("G4ParameterisedNavigation::ComputeStep()", 
                              "GeomNav1002", JustWarning, message);
                  G4cout.precision(oldcoutPrec);
                }
              }
#endif
          }
        }
      }
    }

    if ( initialNode )
    {
      initialNode = false;
      voxelSafety = ComputeVoxelSafety(localPoint,axis);
      if ( voxelSafety<ourSafety )
      {
        ourSafety = voxelSafety;
      }
      if ( currentProposedStepLength<ourSafety )
      {
        // Guaranteed physics limited
        //      
        noStep = false;
        entering = false;
        exiting = false;
        *pBlockedPhysical = nullptr;
        ourStep = kInfinity;
      }
      else
      {
        // Consider intersection with mother solid
        //
        if ( motherSafety<=ourStep )
        {
          if ( !fCheck )           
          {
            motherStep = motherSolid->DistanceToOut(localPoint,
                                                   localDirection,
                                                   true,
                                                   &motherValidExitNormal,
                                                   &motherExitNormal);
          }

          if( ( motherStep < 0.0 ) || ( motherStep >= kInfinity) )
          {
#ifdef G4VERBOSE
            fLogger->ReportOutsideMother(localPoint, localDirection,
                                         motherPhysical);
#endif
            ourStep = motherStep = 0.0;
            // Rely on the code below to set the remaining state, i.e.
            // exiting, entering,  exitNormal & validExitNormal,
            // pBlockedPhysical etc.
          }
#ifdef G4VERBOSE
          if( motherValidExitNormal && ( fCheck || (motherStep<=ourStep)) )
          {
            fLogger->CheckAndReportBadNormal(motherExitNormal,
                                             localPoint, localDirection,
                                             motherStep, motherSolid,
                                             "From motherSolid::DistanceToOut");
          }
#endif
          if ( motherStep<=ourStep )
          {
            ourStep = motherStep;
            exiting = true;
            entering = false;
            if ( validExitNormal )
            {
              const G4RotationMatrix* rot = motherPhysical->GetRotation();
              if (rot)
              {
                exitNormal *= rot->inverse();
              }
            }
          }
          else
          {
            validExitNormal = false;
          }
        }
      }
      newSafety = ourSafety;
    }
    if (noStep)
    {
      noStep = LocateNextVoxel(localPoint, localDirection, ourStep, axis);
    }
  } while (noStep);

  return ourStep;
}

// ***************************************************************************
// ComputeSafety
// ***************************************************************************
//
G4double
G4ParameterisedNavigation::ComputeSafety(const G4ThreeVector& localPoint,
                                         const G4NavigationHistory& history,
                                         const G4double )
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical;
  G4VPVParameterisation *sampleParam;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid, *sampleSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo, curVoxelNodeNo;

  G4SmartVoxelNode *curVoxelNode;
  G4long curNoVolumes, contentNo;
  G4double voxelSafety;

  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;

  motherPhysical = history.GetTopVolume();
  motherLogical = motherPhysical->GetLogicalVolume();
  motherSolid = motherLogical->GetSolid();

  //
  // Compute mother safety
  //

  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety;                     // Working isotropic safety

  //
  // Compute daughter safeties
  //

  // By definition, parameterised volumes exist as first
  // daughter of the mother volume
  //
  samplePhysical = motherLogical->GetDaughter(0);
  samplePhysical->GetReplicationData(axis, nReplicas,
                                     width, offset, consuming);
  sampleParam = samplePhysical->GetParameterisation();

  // Look inside the current Voxel only at the current point
  //
  if ( axis==kUndefined )      // 3D case: current voxel node is retrieved
  {                            //          from G4VoxelNavigation.
    curVoxelNode = fVoxelNode;
  }
  else                         // 1D case: current voxel node is computed here.
  {
    curVoxelNodeNo = G4int((localPoint(fVoxelAxis)
                           -fVoxelHeader->GetMinExtent()) / fVoxelSliceWidth );
    curVoxelNode = fVoxelHeader->GetSlice(curVoxelNodeNo)->GetNode();
    fVoxelNodeNo = curVoxelNodeNo;
    fVoxelNode = curVoxelNode;
  }
  curNoVolumes = curVoxelNode->GetNoContained();

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    sampleNo = curVoxelNode->GetVolume((G4int)contentNo);
    
    // Call virtual methods, and copy information if needed
    //
    sampleSolid= IdentifyAndPlaceSolid( sampleNo,samplePhysical,sampleParam ); 

    G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                               samplePhysical->GetTranslation());
    sampleTf.Invert();
    const G4ThreeVector samplePoint = sampleTf.TransformPoint(localPoint);
    G4double sampleSafety = sampleSolid->DistanceToIn(samplePoint);
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }
  }

  voxelSafety = ComputeVoxelSafety(localPoint,axis);
  if ( voxelSafety<ourSafety )
  {
    ourSafety=voxelSafety;
  }

  return ourSafety;
}

// ********************************************************************
// ComputeVoxelSafety
//
// Computes safety from specified point to collected voxel boundaries
// using already located point.
// ********************************************************************
//
G4double G4ParameterisedNavigation::
ComputeVoxelSafety(const G4ThreeVector& localPoint,
                   const EAxis pAxis) const
{
  // If no best axis is specified, adopt default
  // strategy as for placements
  //  
  if ( pAxis==kUndefined )
  {
    return G4VoxelNavigation::ComputeVoxelSafety(localPoint);
  }

  G4double voxelSafety, plusVoxelSafety, minusVoxelSafety;
  G4double curNodeOffset, minCurCommonDelta, maxCurCommonDelta;
  G4long minCurNodeNoDelta, maxCurNodeNoDelta;
  
  // Compute linear intersection distance to boundaries of max/min
  // to collected nodes at current level
  //
  curNodeOffset = fVoxelNodeNo*fVoxelSliceWidth;
  minCurCommonDelta = localPoint(fVoxelAxis)
                    - fVoxelHeader->GetMinExtent()-curNodeOffset;
  maxCurNodeNoDelta = fVoxelNode->GetMaxEquivalentSliceNo()-fVoxelNodeNo;
  minCurNodeNoDelta = fVoxelNodeNo-fVoxelNode->GetMinEquivalentSliceNo();
  maxCurCommonDelta = fVoxelSliceWidth-minCurCommonDelta;
  plusVoxelSafety   = minCurNodeNoDelta*fVoxelSliceWidth+minCurCommonDelta;
  minusVoxelSafety  = maxCurNodeNoDelta*fVoxelSliceWidth+maxCurCommonDelta;
  voxelSafety = std::min(plusVoxelSafety,minusVoxelSafety);

  if ( voxelSafety<0 )
  {
    voxelSafety = 0;
  }

  return voxelSafety;
}

// ********************************************************************
// LocateNextVoxel
//
// Finds the next voxel from the current voxel and point
// in the specified direction.
//
// Returns false if all voxels considered
//         true  otherwise
// [current Step ends inside same voxel or leaves all voxels]
// ********************************************************************
//
G4bool G4ParameterisedNavigation::
LocateNextVoxel( const G4ThreeVector& localPoint,
                 const G4ThreeVector& localDirection,
                 const G4double currentStep,
                 const EAxis pAxis)
{
  // If no best axis is specified, adopt default
  // location strategy as for placements
  //  
  if ( pAxis==kUndefined )
  {
    return G4VoxelNavigation::LocateNextVoxel(localPoint,
                                              localDirection,
                                              currentStep);
  }

  G4bool isNewVoxel;
  G4int newNodeNo;
  G4double minVal, maxVal, curMinExtent, curCoord;

  curMinExtent = fVoxelHeader->GetMinExtent();
  curCoord = localPoint(fVoxelAxis)+currentStep*localDirection(fVoxelAxis);
  minVal = curMinExtent+fVoxelNode->GetMinEquivalentSliceNo()*fVoxelSliceWidth;
  isNewVoxel = false;

  if ( minVal<=curCoord )
  {
    maxVal = curMinExtent
           + (fVoxelNode->GetMaxEquivalentSliceNo()+1)*fVoxelSliceWidth;
    if ( maxVal<curCoord )
    {
      newNodeNo = fVoxelNode->GetMaxEquivalentSliceNo()+1;
      if ( newNodeNo<G4int(fVoxelHeader->GetNoSlices()) )
      {
        fVoxelNodeNo = newNodeNo;
        fVoxelNode = fVoxelHeader->GetSlice(newNodeNo)->GetNode();
        isNewVoxel = true;
      }
    }
  }
  else
  {
    newNodeNo = fVoxelNode->GetMinEquivalentSliceNo()-1;

    // Must locate from newNodeNo no and down to setup stack and fVoxelNode
    // Repeat or earlier code...
    //
    if ( newNodeNo>=0 )
    {
      fVoxelNodeNo = newNodeNo;
      fVoxelNode = fVoxelHeader->GetSlice(newNodeNo)->GetNode();
      isNewVoxel = true;
    }
  }
  return isNewVoxel;
}

// ********************************************************************
// LevelLocate
// ********************************************************************
//
G4bool
G4ParameterisedNavigation::LevelLocate( G4NavigationHistory& history,
                                  const G4VPhysicalVolume* blockedVol,
                                  const G4int blockedNum,
                                  const G4ThreeVector& globalPoint,
                                  const G4ThreeVector* globalDirection,
                                  const G4bool pLocatedOnEdge, 
                                        G4ThreeVector& localPoint )
{
  G4SmartVoxelHeader *motherVoxelHeader;
  G4SmartVoxelNode *motherVoxelNode;
  G4VPhysicalVolume *motherPhysical, *pPhysical;
  G4VPVParameterisation *pParam;
  G4LogicalVolume *motherLogical;
  G4VSolid *pSolid;
  G4ThreeVector samplePoint;
  G4int voxelNoDaughters, replicaNo;
  
  motherPhysical = history.GetTopVolume();
  motherLogical = motherPhysical->GetLogicalVolume();
  motherVoxelHeader = motherLogical->GetVoxelHeader();

  // Find the voxel containing the point
  //
  motherVoxelNode = ParamVoxelLocate(motherVoxelHeader,localPoint);
  
  voxelNoDaughters = (G4int)motherVoxelNode->GetNoContained();
  if ( voxelNoDaughters==0 )  { return false; }
  
  pPhysical = motherLogical->GetDaughter(0);
  pParam = pPhysical->GetParameterisation();

  // Save parent history in touchable history
  //   ... for use as parent t-h in ComputeMaterial method of param
  //
  G4TouchableHistory parentTouchable( history ); 

  // Search replicated daughter volume
  //
  for ( auto sampleNo=voxelNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    replicaNo = motherVoxelNode->GetVolume(sampleNo);
    if ( (replicaNo!=blockedNum) || (pPhysical!=blockedVol) )
    {
      // Obtain solid (as it can vary) and obtain its parameters
      //
      pSolid = IdentifyAndPlaceSolid( replicaNo, pPhysical, pParam ); 

      // Setup history
      //
      history.NewLevel(pPhysical, kParameterised, replicaNo);
      samplePoint = history.GetTopTransform().TransformPoint(globalPoint);
      if ( !G4AuxiliaryNavServices::CheckPointOnSurface( pSolid, samplePoint,
            globalDirection, history.GetTopTransform(), pLocatedOnEdge) )
      {
        history.BackLevel();
      }
      else
      { 
        // Enter this daughter
        //
        localPoint = samplePoint;
        
        // Set the correct copy number in physical
        //
        pPhysical->SetCopyNo(replicaNo);
        
        // Set the correct solid and material in Logical Volume
        //
        G4LogicalVolume *pLogical = pPhysical->GetLogicalVolume();
        pLogical->SetSolid(pSolid);
        pLogical->UpdateMaterial(pParam->ComputeMaterial(replicaNo,
                                 pPhysical, &parentTouchable)  );
        return true;
      }
    }
  }
  return false;
}
