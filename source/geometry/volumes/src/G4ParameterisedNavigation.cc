//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParameterisedNavigation.cc,v 1.6 2002-05-15 10:23:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4ParameterisedNavigation Implementation
//
// ********************************************************************

#include "G4ParameterisedNavigation.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4ParameterisedNavigation::G4ParameterisedNavigation()
  : fVoxelHeader(0)
{
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4ParameterisedNavigation::~G4ParameterisedNavigation()
{
#ifdef G4DEBUG_NAVIGATION
  G4cout << "G4ParameterisedNavigation::~G4ParameterisedNavigation() called."
   << G4endl;
#endif
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
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int sampleNo;

  G4bool initialNode, noStep;
  G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
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

  //
  // Compute daughter safeties & intersections
  //

  initialNode = true;
  noStep = true;

  // By definition, parameterised volumes exist as first
  // daughter of the mother volume
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
      assert( (0 <= blockedReplicaNo)&&(blockedReplicaNo<nReplicas) );
      //
      // Block exited daughter replica; Must be on boundary => zero safety
      //
      fBList.BlockVolume(blockedReplicaNo);
      ourSafety = 0;
    }
  }
  exiting = false;
  entering = false;

  sampleParam = samplePhysical->GetParameterisation();

  do
  {
    curVoxelNode = fVoxelNode;
    curNoVolumes = curVoxelNode->GetNoContained();

    for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
    {
      sampleNo = curVoxelNode->GetVolume(contentNo);
      if ( !fBList.IsBlocked(sampleNo) )
      {
        fBList.BlockVolume(sampleNo);
        sampleSolid = sampleParam->ComputeSolid(sampleNo, samplePhysical);
        sampleSolid->ComputeDimensions(sampleParam, sampleNo, samplePhysical);
        sampleParam->ComputeTransformation(sampleNo, samplePhysical);
        samplePhysical->Setup(motherPhysical);
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
        *pBlockedPhysical = 0;
        ourStep = kInfinity;
      }
      else
      {
        // Compute mother intersection if required
        //
        if ( motherSafety<=ourStep )
        {
          G4double motherStep = motherSolid->DistanceToOut(localPoint,
                                                           localDirection,
                                                           true,
                                                           &validExitNormal,
                                                           &exitNormal);
          if ( motherStep<=ourStep )
          {
            ourStep = motherStep;
            exiting = true;
            entering = false;
            if ( validExitNormal )
            {
              const G4RotationMatrix *rot = motherPhysical->GetRotation();
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
      newSafety=ourSafety;
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
                                         const G4double pProposedMaxLength )
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical;
  G4VPVParameterisation *sampleParam;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid, *sampleSolid;
  G4double motherSafety, ourSafety;
  G4int sampleNo, curVoxelNodeNo;

  G4SmartVoxelNode *curVoxelNode;
  G4int curNoVolumes, contentNo;
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
  samplePhysical->GetReplicationData(axis, nReplicas, width, offset, consuming);
  sampleParam = samplePhysical->GetParameterisation();

  // Look inside the current Voxel only at the current point
  //
  if ( axis==kUndefined )      // 3D case: current voxel node is retrieved
  {                            //          from G4VoxelNavigation.
    curVoxelNode = fVoxelNode;
  }
  else                         // 1D case: current voxel node is computed here.
  {
    curVoxelNodeNo = G4int((localPoint(fVoxelAxis)-fVoxelHeader->GetMinExtent())
                           / fVoxelSliceWidth );
    curVoxelNode = fVoxelHeader->GetSlice(curVoxelNodeNo)->GetNode();
    fVoxelNodeNo = curVoxelNodeNo;
    fVoxelNode = curVoxelNode;
  }
  curNoVolumes = curVoxelNode->GetNoContained();

  for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
  {
    sampleNo = curVoxelNode->GetVolume(contentNo);
    sampleSolid = sampleParam->ComputeSolid(sampleNo, samplePhysical);
    sampleSolid->ComputeDimensions(sampleParam, sampleNo, samplePhysical);
    sampleParam->ComputeTransformation(sampleNo, samplePhysical);
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
    return G4VoxelNavigation::ComputeVoxelSafety(localPoint);

  G4double voxelSafety, plusVoxelSafety, minusVoxelSafety;
  G4double curNodeOffset, minCurCommonDelta, maxCurCommonDelta;
  G4int minCurNodeNoDelta, maxCurNodeNoDelta;
  
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
  voxelSafety = G4std::min(plusVoxelSafety,minusVoxelSafety);

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
    return G4VoxelNavigation::LocateNextVoxel(localPoint,
                                              localDirection,
                                              currentStep);
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
      if ( newNodeNo<fVoxelHeader->GetNoSlices() )
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
