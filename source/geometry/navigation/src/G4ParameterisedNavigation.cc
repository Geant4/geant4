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
// $Id: G4ParameterisedNavigation.cc,v 1.13 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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

// ********************************************************************
// Constructor
// ********************************************************************
//
G4ParameterisedNavigation::G4ParameterisedNavigation()
  : fVoxelAxis(kUndefined), fVoxelNoSlices(0), fVoxelSliceWidth(0.),
    fVoxelNodeNo(0), fVoxelHeader(0)
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

#ifdef G4VERBOSE
  if ( fCheck )
  {
    if( motherSafety < 0.0 )
    {
      G4cout << "ERROR - G4ParameterisedNavigation::ComputeStep()" << G4endl
             << "        Current solid " << motherSolid->GetName()
             << " gave negative safety: " << motherSafety << G4endl
             << "        for the current (local) point " << localPoint
             << G4endl;
      motherSolid->DumpInfo();
      G4Exception("G4ParameterisedNavigation::ComputeStep()",
                  "NegativeSafetyMotherVol", FatalException,
                  "Negative Safety In Voxel Navigation !" ); 
    }
    if( motherSolid->Inside(localPoint)==kOutside )
    { 
      G4cout << "WARNING - G4ParameterisedNavigation::ComputeStep()" << G4endl
             << "          Point " << localPoint
             << " is outside current volume " << motherPhysical->GetName()
             << G4endl;
      G4double  estDistToSolid= motherSolid->DistanceToIn(localPoint); 
      G4cout << "          Estimated isotropic distance to solid (distToIn)= " 
             << estDistToSolid << G4endl;
      if( estDistToSolid > 100.0 * motherSolid->GetTolerance() )
      {
        motherSolid->DumpInfo();
        G4Exception("G4ParameterisedNavigation::ComputeStep()",
                    "FarOutsideCurrentVolume", FatalException,
                    "Point is far outside Current Volume !"); 
      }
      else
        G4Exception("G4ParameterisedNavigation::ComputeStep()",
                    "OutsideCurrentVolume", JustWarning,
       "Point is a little outside Current Volume."); 
    }
  }
#endif

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
                intersectionPoint= samplePoint + sampleStep * sampleDirection;
                EInside insideIntPt= sampleSolid->Inside(intersectionPoint); 
                if( insideIntPt != kSurface )
                {
                  G4int oldcoutPrec = G4cout.precision(16); 
                  G4cout << "WARNING - G4ParameterisedNavigation::ComputeStep()"
                         << G4endl
                         << "          Inaccurate solid DistanceToIn"
                         << " for solid " << sampleSolid->GetName() << G4endl;
                  G4cout << "          Solid gave DistanceToIn = "
                         << sampleStep << " yet returns " ;
                  if( insideIntPt == kInside )
                    G4cout << "-kInside-"; 
                  else if( insideIntPt == kOutside )
                    G4cout << "-kOutside-";
                  else
                    G4cout << "-kSurface-"; 
                  G4cout << " for this point !" << G4endl; 
                  G4cout << "          Point = " << intersectionPoint << G4endl;
                  if ( insideIntPt != kInside )
                    G4cout << "        DistanceToIn(p) = " 
                           << sampleSolid->DistanceToIn(intersectionPoint)
                           << G4endl;
                  if ( insideIntPt != kOutside ) 
                    G4cout << "        DistanceToOut(p) = " 
                           << sampleSolid->DistanceToOut(intersectionPoint)
                           << G4endl;
                  G4Exception("G4ParameterisedNavigation::ComputeStep()", 
                              "InaccurateDistanceToIn", JustWarning,
                              "Navigator gets conflicting response from Solid.");
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
        *pBlockedPhysical = 0;
        ourStep = kInfinity;
      }
      else
      {
        //
        // Compute mother intersection if required
        //
        if ( motherSafety<=ourStep )
        {
          G4double motherStep = motherSolid->DistanceToOut(localPoint,
                                                           localDirection,
                                                           true,
                                                           &validExitNormal,
                                                           &exitNormal);
#ifdef G4VERBOSE
          if ( fCheck ) 
            if( ( motherStep < 0.0 ) || ( motherStep >= kInfinity) )
            {
              G4int oldPrOut= G4cout.precision(16); 
              G4int oldPrErr= G4cerr.precision(16);
              G4cerr << "ERROR - G4ParameterisedNavigation::ComputeStep()"
                     << G4endl
                     << "        Problem in Navigation"  << G4endl
                     << "        Point (local coordinates): "
                     << localPoint << G4endl
                     << "        Local Direction: "
                     << localDirection << G4endl
                     << "        Solid: " << motherSolid->GetName() << G4endl; 
              motherSolid->DumpInfo();
              G4Exception("G4ParameterisedNavigation::ComputeStep()",
                          "PointOutsideCurrentVolume", FatalException,
                          "Current point is outside the current solid !");
              G4cout.precision(oldPrOut);
              G4cerr.precision(oldPrErr);
            }
#endif
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
                                         const G4double )
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
    sampleNo = curVoxelNode->GetVolume(contentNo);
    
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
  
  voxelNoDaughters = motherVoxelNode->GetNoContained();
  if ( voxelNoDaughters==0 )  { return false; }
  
  pPhysical = motherLogical->GetDaughter(0);
  pParam = pPhysical->GetParameterisation();

  // Save parent history in touchable history
  //   ... for use as parent t-h in ComputeMaterial method of param
  //
  G4TouchableHistory parentTouchable( history ); 

  // Search replicated daughter volume
  //
  for ( register int sampleNo=voxelNoDaughters-1; sampleNo>=0; sampleNo-- )
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
