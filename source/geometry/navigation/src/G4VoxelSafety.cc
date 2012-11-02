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
// $Id: G4VoxelSafety.cc,v 1.9 2010-11-11 16:15:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Author:  John Apostolakis
//  First version:  31 May 2010
// 
// --------------------------------------------------------------------
#include "G4VoxelSafety.hh"

#include "G4GeometryTolerance.hh"

#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"

// State 
//    - values used during computation of Safety 
// 
// Constructor
//     - copied from G4VoxelNavigation  (1st version)
G4VoxelSafety::G4VoxelSafety()
  : fBlockList(),
    fpMotherLogical(0),
    fVoxelDepth(-1),
    fVoxelAxisStack(kNavigatorVoxelStackMax,kXAxis),
    fVoxelNoSlicesStack(kNavigatorVoxelStackMax,0),
    fVoxelSliceWidthStack(kNavigatorVoxelStackMax,0.),
    fVoxelNodeNoStack(kNavigatorVoxelStackMax,0),
    fVoxelHeaderStack(kNavigatorVoxelStackMax,(G4SmartVoxelHeader*)0),
    fVoxelNode(0),
    fCheck(false),
    fVerbose(0)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
   // fVerbose=3;
}

G4VoxelSafety::~G4VoxelSafety()
{
}

// ********************************************************************
// ComputeSafety
//
// Calculates the isotropic distance to the nearest boundary from the
// specified point in the local coordinate system. 
// The localpoint utilised must be within the current volume.
// ********************************************************************
//
G4double
G4VoxelSafety::ComputeSafety(const G4ThreeVector&     localPoint,
                             const G4VPhysicalVolume& currentPhysical,
                                   G4double ) //          maxLength)
{
  // G4VPhysicalVolume *samplePhysical;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4SmartVoxelHeader *motherVoxelHeader;
  G4double motherSafety, ourSafety;
  G4int localNoDaughters;  // , sampleNo;
  // G4SmartVoxelNode *curVoxelNode;
  // G4int    curNoVolumes, contentNo;
  G4double daughterSafety;

  motherLogical = currentPhysical.GetLogicalVolume();
  fpMotherLogical= motherLogical;   // For use by the other methods
  motherSolid = motherLogical->GetSolid();
  motherVoxelHeader= motherLogical->GetVoxelHeader();

#ifdef G4VERBOSE  
  if( fVerbose > 0 )
  { 
    G4cout << "*** G4VoxelSafety::ComputeSafety(): ***" << G4endl; 
  }
#endif

  // Check that point is inside mother volume
  EInside  insideMother= motherSolid->Inside(localPoint); 
  if( insideMother != kInside  )
  { 
#ifdef G4DEBUG_NAVIGATION
    if( insideMother == kOutside )
    {
      std::ostringstream message;
      message << "Safety method called for location outside current Volume." << G4endl
              << "Location for safety is Outside this volume. " << G4endl
              << "The approximate distance to the solid "
              << "(safety from outside) is: " 
              << motherSolid->DistanceToIn( localPoint ) << G4endl;
      message << "  Problem occurred with physical volume: "
              << " Name: " << currentPhysical.GetName()
              << " Copy No: " << currentPhysical.GetCopyNo() << G4endl
              << "    Local Point = " << localPoint << G4endl;
      message << "  Description of solid: " << G4endl
              << *motherSolid << G4endl;
      G4Exception("G4VoxelSafety::ComputeSafety()", "GeomNav0003",
                  FatalException, message);
    }
#endif
    return 0.0;
  }   

  //  First limit:  mother safety - distance to outer boundaries
  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety
   
  // fCheck= 1;  //  <- Use this to check / debug only 
#ifdef G4VERBOSE
  if(( fCheck ) ) // && ( fVerbose == 1 ))
  {
    // G4cout << "*** G4VoxelSafety::ComputeSafety(): ***" << G4endl;
    G4cout << "    Invoked DistanceToOut(p) for mother solid: "
           << motherSolid->GetName()
           << ". Solid replied: " << motherSafety << G4endl
           << "    For local point p: " << localPoint
           << ", to be considered as 'mother safety'." << G4endl;
  }
#endif
  localNoDaughters = motherLogical->GetNoDaughters();

  fBlockList.Enlarge(localNoDaughters);
  fBlockList.Reset();

  fVoxelDepth = -1;  // Resets the depth -- must be done for next method
  daughterSafety= SafetyForVoxelHeader( motherVoxelHeader, localPoint ); 

  ourSafety= std::min( motherSafety, daughterSafety ); 

  return ourSafety;
}

// Calculate the safety for volumes included in current Voxel Node
// 
G4double
G4VoxelSafety::
SafetyForVoxelNode( const G4SmartVoxelNode *curVoxelNode,
                    const G4ThreeVector&   localPoint )
{
   G4double ourSafety= DBL_MAX;

   G4int    curNoVolumes, contentNo, sampleNo;
   G4VPhysicalVolume *samplePhysical;

   G4double      sampleSafety=0.0; 
   G4ThreeVector samplePoint;
   G4VSolid*     ptrSolid=0;

   curNoVolumes = curVoxelNode->GetNoContained();

   for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
   {
      sampleNo = curVoxelNode->GetVolume(contentNo);
      if ( !fBlockList.IsBlocked(sampleNo) ) 
      { 
        fBlockList.BlockVolume(sampleNo);

        samplePhysical = fpMotherLogical->GetDaughter(sampleNo);
        G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                   samplePhysical->GetTranslation());
        sampleTf.Invert();
        samplePoint =  sampleTf.TransformPoint(localPoint);
        ptrSolid =  samplePhysical->GetLogicalVolume()->GetSolid();

        sampleSafety = ptrSolid->DistanceToIn(samplePoint);
        ourSafety   = std::min( sampleSafety, ourSafety ); 
#ifdef G4VERBOSE
        if(( fCheck ) && ( fVerbose == 1 ))
        {
          // ReportSolidSafetyToIn( MethodName, solid, value, point ); 
          G4cout << "*** G4VoxelSafety::SafetyForVoxelNode(): ***" << G4endl
                 << "    Invoked DistanceToIn(p) for daughter solid: "
                 << ptrSolid->GetName()
                 << ". Solid replied: " << sampleSafety << G4endl
                 << "    For local point p: " << samplePoint
                 << ", to be considered as 'daughter safety'." << G4endl;
        }
#endif
      }
   }  // end for contentNo

   return ourSafety; 
}


//
//  Pseudo-code for  Compute Safety
//
//  for each (potential) dimension of depth
//     iterate through VoxelProxies (Header, Node) 
//      while  distanceToVoxel  <= estimatedSafety
// 
//  Better:
//    iterate through/down the tree of VoxelProxies 
//      while    distanceToVoxel <= estimatedSafety
//      Examine each node at current level for which this condition holds. 
//                   (version 0 can examine all nodes)

// ********************************************************************
//  SafetyForVoxelHeader method
//   - which cycles through levels of headers to process each node level
//   - Obtained by modifying VoxelLocate (to cycle through Node Headers)
// *********************************************************************
//
G4double
G4VoxelSafety::SafetyForVoxelHeader( const G4SmartVoxelHeader* pHeader,
                                     const G4ThreeVector&   localPoint,
                                           G4double     distUpperDepth )
{
  const G4SmartVoxelHeader * const targetVoxelHeader=pHeader;
  G4SmartVoxelNode *targetVoxelNode=0;

  const G4SmartVoxelProxy *sampleProxy;
  EAxis    targetHeaderAxis;
  G4double targetHeaderMin, targetHeaderMax, targetHeaderNodeWidth;
  G4int    targetHeaderNoSlices;
  G4int    targetNodeNo;   // ,  pointNodeNo;
  // G4int    minCurNodeNoDelta, maxCurNodeNoDelta;
  static const char *pMethodName= "G4VoxelSafety::SafetyForVoxelHeader";
  G4double minSafety= DBL_MAX;
  unsigned int checkedNum= 0;
  
  fVoxelDepth++;
  // fVoxelDepth  set by ComputeSafety or previous level call

  targetHeaderAxis =      targetVoxelHeader->GetAxis();
  targetHeaderNoSlices =  targetVoxelHeader->GetNoSlices();
  targetHeaderMin =       targetVoxelHeader->GetMinExtent();
  targetHeaderMax =       targetVoxelHeader->GetMaxExtent();

  targetHeaderNodeWidth = (targetHeaderMax-targetHeaderMin)
                          / targetHeaderNoSlices;

  G4double localCrd= localPoint(targetHeaderAxis);

  const G4int pointNodeNo= G4int( (localCrd-targetHeaderMin)
                                / targetHeaderNodeWidth );
  // Ensure that it is between 0 and targetHeader->GetMaxExtent() - 1

#ifdef G4VERBOSE  
  if( fVerbose > 2 )
  { 
    G4cout << G4endl;
    G4cout << "**** G4VoxelSafety::SafetyForVoxelHeader  " << G4endl;
    G4cout << "  Called at Depth     = " << fVoxelDepth ; // << G4endl;
    G4cout << " Calculated pointNodeNo= " << pointNodeNo
           << "  from position= " <<  localPoint(targetHeaderAxis)
           << "  min= "    << targetHeaderMin
           << "  max= "    << targetVoxelHeader->GetMaxExtent()
           << "  width= "  << targetHeaderNodeWidth 
           << "  no-slices= " << targetHeaderNoSlices
           << "  axis=  "  << targetHeaderAxis   << G4endl;
  }else if (fVerbose == 1){
    G4cout << " VoxelSafety: Depth  = " << fVoxelDepth 
           << " Number of Slices = " << targetHeaderNoSlices
           << " Header (address) = " << targetVoxelHeader  << G4endl;
  }
#endif

  // Stack info for stepping
  //
  fVoxelAxisStack[fVoxelDepth] = targetHeaderAxis;
  fVoxelNoSlicesStack[fVoxelDepth] = targetHeaderNoSlices;
  fVoxelSliceWidthStack[fVoxelDepth] = targetHeaderNodeWidth;

  fVoxelHeaderStack[fVoxelDepth] = pHeader;


  // G4int  numSlicesCheck= targetHeaderNoSlices; 

  // G4cout << "---> Current Voxel Header has " << *targetVoxelHeader << G4endl;

  G4int nextUp,   trialUp= -1;    // =   pointNodeNo+1;
  G4int nextDown, trialDown= -1;  // = pointNodeNo-1;

  G4double distUp= DBL_MAX, distDown= DBL_MAX;

   
// #ifdef SIMPLE
  // Ignore equivalents for now
  nextUp=   pointNodeNo+1;
  nextDown= pointNodeNo-1;

  // Use Equivalent voxels -- NO this should not be at this level - it should be at level below!!!
  // nextUp =   targetVoxelHeader->GetMaxEquivalentSliceNo()+1;
  // nextDown = targetVoxelHeader->GetMinEquivalentSliceNo()-1;
  
  G4int    nextNodeNo= pointNodeNo;
  G4double distAxis;  // Distance in current Axis
  distAxis= 0.0;  // Starting in node containing local Coordinate

  //  for( targetNodeNo= pointNodeNo;
  //               //   (targetNodeNo<targetHeaderNoSlices) &&
  //     (targetNodeNo>=0); targetNodeNo= nextNodeNo )
  G4bool nextIsInside= false;
  targetNodeNo= pointNodeNo;
  do
  {
     G4double nodeSafety= DBL_MAX, levelSafety= DBL_MAX;
     fVoxelNodeNoStack[fVoxelDepth] = targetNodeNo;
    
     checkedNum++; 
    
     sampleProxy = targetVoxelHeader->GetSlice(targetNodeNo);

#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 2 ) {
       G4cout << " -Checking node " << targetNodeNo
              << " is proxy with address " << sampleProxy << G4endl;
     }
#endif 

     if ( sampleProxy == 0 )
     {
       G4ExceptionDescription ed;
       ed << " Problem for node number= " << targetNodeNo
          << "    Number of slides = " << targetHeaderNoSlices
          << G4endl;
       G4Exception( pMethodName,  "GeomNav003", FatalException,
                    ed, "Problem sampleProxy is Zero. Failure in loop.");
       exit(1);
     }
     else if ( sampleProxy->IsNode() )
     {
        targetVoxelNode = sampleProxy->GetNode();

        // Deal with the node here [ i.e. the last level ] 
        nodeSafety= SafetyForVoxelNode( targetVoxelNode, localPoint);
#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2 ) {
           G4cout << " -- It is a Node "; //  << G4endl;
           G4cout << "    its safety= " << nodeSafety
            << "  Min= " << minSafety << G4endl;
        }
#endif
        minSafety= std::min( minSafety, nodeSafety );

        trialUp =   targetVoxelNode->GetMaxEquivalentSliceNo()+1;
        trialDown = targetVoxelNode->GetMinEquivalentSliceNo()-1;
     }
     else  
     {
        const G4SmartVoxelHeader *pNewVoxelHeader = sampleProxy->GetHeader();
        // fVoxelDepth++;

        G4double distCombined;
        distCombined = std::max( distUpperDepth, distAxis);
        // distCombined = std::sqrt( distUpperDepth*distUpperDepth + distAxis*distAxis);
#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2 ) {
          G4cout << " -- It is a Header " << G4endl;
          G4cout << "  Recurse to deal with next level, fVoxelDepth= " 
                 << fVoxelDepth+1 << G4endl;
          G4cout << "  Distances:  Upper= " << distUpperDepth
                 << " Axis= " << distAxis
                 << " Combined= " << distCombined << G4endl;
        }
#endif
        
        // Recurse to deal with lower levels
        levelSafety= SafetyForVoxelHeader( pNewVoxelHeader, localPoint, distCombined);
        // fVoxelDepth--;

#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2 ) { 
          G4cout << "      >> Level safety = " << levelSafety << G4endl;
        }
#endif 
        minSafety= std::min( minSafety, levelSafety );
        trialUp =   pNewVoxelHeader->GetMaxEquivalentSliceNo()+1;
        trialDown = pNewVoxelHeader->GetMinEquivalentSliceNo()-1;
     }

     // Find next closest Voxel 
     //    - first try: by simple subtraction
     //    - later:  using distance  (TODO - tbc)

     if( targetNodeNo >= pointNodeNo )
     {
         nextUp = trialUp;
         distUp = targetHeaderMax-localCrd;
#ifdef G4DEBUG_NAVIGATION
         if( fVerbose > 2 ) {
            G4cout << " > Updated nextUp= " << nextUp << G4endl;
         }
#endif
     }
     
     if( targetNodeNo <= pointNodeNo ) 
     {
         nextDown = trialDown;
         distDown = localCrd-targetHeaderMin;
#ifdef G4DEBUG_NAVIGATION
         if( fVerbose > 2 ) {
            G4cout << " > Updated nextDown= " << nextDown << G4endl;
         }
#endif
     }

#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 2 )   {
        G4cout << " Node= " << pointNodeNo
               << " Up:   next= " << nextUp  << " d# " << nextUp - pointNodeNo
               << " trialUp=  " << trialUp << " d# " << trialUp - pointNodeNo   // << G4endl
               << " Down: next= " << nextDown << " d# " << targetNodeNo - nextDown
               << " trialDwn= " << trialDown << " d# " << targetNodeNo - trialDown // << G4endl
               << " condition (next is Inside)= " << nextIsInside
               << G4endl;
     }
#endif     
     
     G4bool UpIsClosest;
     // UpIsClosest= (nextUp - pointNodeNo) < (pointNodeNo - nextDown);
     UpIsClosest= distUp < distDown;
     // G4cout << " Distances:  Up " << distUp << "  Down= " << distDown
     //        << " Depth = " << fVoxelDepth << G4endl;
     
     if( UpIsClosest || (nextDown < 0) )
     {
        nextNodeNo=nextUp;
        distAxis = distUp;
        ++nextUp; // Default
                  // #ifdef G4VERBOSE
        if( fVerbose > 2 ) {
           G4cout << " > Chose Up.    Nodes: next= " << nextNodeNo
                  << " new nextUp=   " << nextUp
                  << " Dist = " << distAxis << G4endl;
        }
        // #endif
     }
     else
     {
        nextNodeNo=nextDown;
        distAxis = distDown;
        --nextDown; // A default value
                    // #ifdef G4VERBOSE
        if( fVerbose > 2 ) {
           G4cout << " > Chose Down.   Nodes: next= " << nextNodeNo
                  << " new nextDown= " << nextDown
                  << " Dist = " << distAxis << G4endl;
        }
        // #endif
     }

     nextIsInside = (nextNodeNo >= 0) && (nextNodeNo < targetHeaderNoSlices);
     if( nextIsInside )
     {
        targetNodeNo= nextNodeNo;
        assert( targetVoxelHeader->GetSlice(targetNodeNo) != 0 ); 

#ifdef G4DEBUG_NAVIGATION
        G4bool bContinue= (distAxis<minSafety);
        if( !bContinue ){
           if( fVerbose > 2 ){
              G4cout << " Can skip remaining at depth " << targetHeaderAxis // << G4endl
                     << " >>  distAxis= " << distAxis << " minSafety= " << minSafety << G4endl;
           }
        }

        const G4SmartVoxelProxy* nextProxy= targetVoxelHeader->GetSlice(targetNodeNo);
        if ( nextProxy == 0 )
        {
          G4ExceptionDescription ed;
          ed << " Problem for next node number= " << targetNodeNo
             << "    Number of slides = " << targetHeaderNoSlices << G4endl;
          G4Exception( pMethodName,  "GeomNav003", FatalException,
                       ed, "Problem nextProxy is Zero. Failure in loop.");
          exit(1);
        }

#endif
     }
     else
     {
#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2) {
           G4cout << " VoxSaf> depth= " << fVoxelDepth << G4endl;
           G4cout << " VoxSaf> No more candidates:  nodeDown= " << nextDown
                  << " nodeUp= " << nextUp
                  << " lastSlice= " << targetHeaderNoSlices << G4endl;
        }
#endif
     }
  } while ( nextIsInside
           && ( distAxis*distAxis + distUpperDepth*distUpperDepth < minSafety*minSafety ) );

#ifdef G4VERBOSE
  if( fVerbose > 0 )
  { 
    G4cout << " Ended for targetNodeNo -- checked " << checkedNum << " out of "
           << targetHeaderNoSlices << " slices." << G4endl;
    G4cout << " ===== Returning from SafetyForVoxelHeader " 
           << "  Depth= " << fVoxelDepth << G4endl
           << G4endl;  
  }
#endif

  // Go back one level
  fVoxelDepth--; 
  
  return minSafety;
}
