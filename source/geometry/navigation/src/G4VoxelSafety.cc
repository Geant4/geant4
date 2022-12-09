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
//  Author:  John Apostolakis
//  First version:  31 May 2010
// 
// --------------------------------------------------------------------
#include "G4VoxelSafety.hh"

#include "G4GeometryTolerance.hh"

#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"

// ********************************************************************
// Constructor
//     - copied from G4VoxelNavigation  (1st version)
// ********************************************************************
//
G4VoxelSafety::G4VoxelSafety()
  : fBlockList(),
    fVoxelAxisStack(kNavigatorVoxelStackMax,kXAxis),
    fVoxelNoSlicesStack(kNavigatorVoxelStackMax,0),
    fVoxelSliceWidthStack(kNavigatorVoxelStackMax,0.),
    fVoxelNodeNoStack(kNavigatorVoxelStackMax,0),
    fVoxelHeaderStack(kNavigatorVoxelStackMax,(G4SmartVoxelHeader*)nullptr)
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

// ********************************************************************
// Destructor
// ********************************************************************
//
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
G4VoxelSafety::ComputeSafety(const G4ThreeVector& localPoint,
                             const G4VPhysicalVolume& currentPhysical,
                                   G4double maxLength)
{
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4SmartVoxelHeader *motherVoxelHeader;
  G4double motherSafety, ourSafety;
  G4int localNoDaughters;
  G4double daughterSafety;

  motherLogical = currentPhysical.GetLogicalVolume();
  fpMotherLogical= motherLogical;   // For use by the other methods
  motherSolid = motherLogical->GetSolid();
  motherVoxelHeader = motherLogical->GetVoxelHeader();

#ifdef G4VERBOSE  
  if( fVerbose > 0 )
  { 
    G4cout << "*** G4VoxelSafety::ComputeSafety(): ***" << G4endl; 
  }
#endif

  // Check that point is inside mother volume
  //
  EInside  insideMother = motherSolid->Inside(localPoint); 
  if( insideMother != kInside  )
  { 
#ifdef G4DEBUG_NAVIGATION
    if( insideMother == kOutside )
    {
      std::ostringstream message;
      message << "Safety method called for location outside current Volume."
              << G4endl
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
  //
  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety;                 // Working isotropic safety
   
#ifdef G4VERBOSE
  if(( fCheck ) ) // && ( fVerbose == 1 ))
  {
    G4cout << "    Invoked DistanceToOut(p) for mother solid: "
           << motherSolid->GetName()
           << ". Solid replied: " << motherSafety << G4endl
           << "    For local point p: " << localPoint
           << ", to be considered as 'mother safety'." << G4endl;
  }
#endif
  localNoDaughters = (G4int)motherLogical->GetNoDaughters();

  fBlockList.Enlarge(localNoDaughters);
  fBlockList.Reset();

  fVoxelDepth = -1;  // Resets the depth -- must be done for next method
  daughterSafety= SafetyForVoxelHeader(motherVoxelHeader, localPoint, maxLength,
                                       currentPhysical, 0.0, ourSafety);
  ourSafety= std::min( motherSafety, daughterSafety ); 

  return ourSafety;
}

// ********************************************************************
// SafetyForVoxelNode
//
// Calculate the safety for volumes included in current Voxel Node
// ********************************************************************
// 
G4double
G4VoxelSafety::SafetyForVoxelNode( const G4SmartVoxelNode* curVoxelNode,
                                   const G4ThreeVector& localPoint )  
{
   G4double ourSafety = DBL_MAX;

   G4long curNoVolumes, contentNo;
   G4int sampleNo;
   G4VPhysicalVolume* samplePhysical;

   G4double      sampleSafety = 0.0; 
   G4ThreeVector samplePoint;
   G4VSolid*     ptrSolid = nullptr;

   curNoVolumes = curVoxelNode->GetNoContained();

   for ( contentNo=curNoVolumes-1; contentNo>=0; contentNo-- )
   {
      sampleNo = curVoxelNode->GetVolume((G4int)contentNo);
      if ( !fBlockList.IsBlocked(sampleNo) ) 
      { 
        fBlockList.BlockVolume(sampleNo);

        samplePhysical = fpMotherLogical->GetDaughter(sampleNo);
        G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                   samplePhysical->GetTranslation());
        sampleTf.Invert();
        samplePoint = sampleTf.TransformPoint(localPoint);
        ptrSolid = samplePhysical->GetLogicalVolume()->GetSolid();

        sampleSafety = ptrSolid->DistanceToIn(samplePoint);
        ourSafety = std::min( sampleSafety, ourSafety ); 
#ifdef G4VERBOSE
        if(( fCheck ) && ( fVerbose == 1 ))
        {
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

// ********************************************************************
// SafetyForVoxelHeader
//
// Cycles through levels of headers to process each node level
// Obtained by modifying VoxelLocate (to cycle through Node Headers)
// *********************************************************************
//
G4double
G4VoxelSafety::SafetyForVoxelHeader( const G4SmartVoxelHeader* pHeader,
                                     const G4ThreeVector& localPoint,
                                           G4double maxLength,
                                     const G4VPhysicalVolume&  currentPhysical, //Debug
                                           G4double distUpperDepth_Sq,
                                           G4double previousMinSafety
                                   )
{
  const G4SmartVoxelHeader* const targetVoxelHeader = pHeader;
  G4SmartVoxelNode* targetVoxelNode = nullptr;

  const G4SmartVoxelProxy* sampleProxy;
  EAxis    targetHeaderAxis;
  G4double targetHeaderMin, targetHeaderMax, targetHeaderNodeWidth;
  G4int    targetHeaderNoSlices;
  G4int    targetNodeNo;

  G4double minSafety = previousMinSafety;
  G4double ourSafety = DBL_MAX;
  unsigned int checkedNum= 0;
  
  ++fVoxelDepth;
  // fVoxelDepth  set by ComputeSafety or previous level call

  targetHeaderAxis =      targetVoxelHeader->GetAxis();
  targetHeaderNoSlices =  (G4int)targetVoxelHeader->GetNoSlices();
  targetHeaderMin =       targetVoxelHeader->GetMinExtent();
  targetHeaderMax =       targetVoxelHeader->GetMaxExtent();

  targetHeaderNodeWidth = (targetHeaderMax-targetHeaderMin)
                          / targetHeaderNoSlices;

  G4double localCrd = localPoint(targetHeaderAxis);

  const G4int candNodeNo = G4int( (localCrd-targetHeaderMin)
                                 / targetHeaderNodeWidth );
  // Ensure that it is between 0 and targetHeader->GetMaxExtent() - 1

#ifdef G4DEBUG_VOXELISATION  
  if( candNodeNo < 0 || candNodeNo > targetHeaderNoSlices-1 )
  {
     G4ExceptionDescription ed;
     ed << " Potential ERROR."
        << " Point is outside range of Voxel in current coordinate" << G4endl;
     ed << "  Node number of point " << localPoint
        << "is outside the range. " << G4endl;
     ed << "    Voxel node Num= " << candNodeNo << " versus  minimum= " << 0
        << " and maximum= " << targetHeaderNoSlices-1 << G4endl;
     ed << "    Axis = " << targetHeaderAxis
        << "  No of slices = " << targetHeaderNoSlices << G4endl;
     ed << "    Local coord = " << localCrd
        << "  Voxel Min = " << targetHeaderMin
        <<            " Max = " << targetHeaderMax << G4endl;
     G4LogicalVolume *pLogical= currentPhysical.GetLogicalVolume();
     ed << "  Current volume (physical) = " << currentPhysical.GetName()
        << "   (logical) = " << pLogical->GetName()     << G4endl;
     G4VSolid* pSolid= pLogical->GetSolid();
     ed << "  Solid type = " << pSolid->GetEntityType() << G4endl;
     ed << *pSolid << G4endl;
     G4Exception("G4VoxelSafety::SafetyForVoxelHeader()", "GeomNav1003",
                 JustWarning, ed,
                 "Point is outside range of Voxel in current coordinate");
  }
#endif
  
  const G4int pointNodeNo =
              std::max( 0, std::min( candNodeNo, targetHeaderNoSlices-1 ) );
  
#ifdef G4VERBOSE  
  if( fVerbose > 2 )
  { 
    G4cout << G4endl;
    G4cout << "**** G4VoxelSafety::SafetyForVoxelHeader  " << G4endl;
    G4cout << "  Called at Depth     = " << fVoxelDepth;
    G4cout << " Calculated pointNodeNo= " << pointNodeNo
           << "  from position= " <<  localPoint(targetHeaderAxis)
           << "  min= "    << targetHeaderMin
           << "  max= "    << targetHeaderMax
           << "  width= "  << targetHeaderNodeWidth 
           << "  no-slices= " << targetHeaderNoSlices
           << "  axis=  "  << targetHeaderAxis   << G4endl;
  }
  else if (fVerbose == 1)
  {
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

  G4int   trialUp = -1,     trialDown = -1;
  G4double distUp = DBL_MAX, distDown = DBL_MAX;

  // Using Equivalent voxels - this is pre-initialisation only
  //
  G4int nextUp =   pointNodeNo+1;
  G4int nextDown = pointNodeNo-1;

  G4int    nextNodeNo = pointNodeNo;
  G4double distAxis;  // Distance in current Axis
  distAxis = 0.0;  // Starting in node containing local Coordinate

  G4bool nextIsInside = false;

  G4double distMaxInterest= std::min( previousMinSafety, maxLength);
    // We will not look beyond this distance.
    // This distance will be updated to reflect the
    // max ( minSafety, maxLength ) at each step

  targetNodeNo = pointNodeNo;
  do
  {
     G4double nodeSafety = DBL_MAX, headerSafety = DBL_MAX;
     fVoxelNodeNoStack[fVoxelDepth] = targetNodeNo;
    
     ++checkedNum; 
    
     sampleProxy = targetVoxelHeader->GetSlice(targetNodeNo);

#ifdef G4DEBUG_NAVIGATION
     assert( sampleProxy != 0);
     if( fVerbose > 2 ) 
     {
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
       G4Exception( "G4VoxelSafety::SafetyForVoxelHeader()", "GeomNav0003",
                    FatalException, ed,
                    "Problem sampleProxy is Zero. Failure in loop.");
     }
     else if ( sampleProxy->IsNode() )
     {
       targetVoxelNode = sampleProxy->GetNode();

       // Deal with the node here [ i.e. the last level ]
       //
       nodeSafety= SafetyForVoxelNode( targetVoxelNode, localPoint);
#ifdef G4DEBUG_NAVIGATION
       if( fVerbose > 2 )
       {
          G4cout << " -- It is a Node ";
          G4cout << "    its safety= " << nodeSafety
                 << " our level Saf = " << ourSafety
                 << "  MinSaf= " << minSafety << G4endl;
       }
#endif
        ourSafety= std::min( ourSafety, nodeSafety );
       
        trialUp = targetVoxelNode->GetMaxEquivalentSliceNo()+1;
        trialDown = targetVoxelNode->GetMinEquivalentSliceNo()-1;
     }
     else  
     {
        const G4SmartVoxelHeader* pNewVoxelHeader = sampleProxy->GetHeader();

        G4double distCombined_Sq;
        distCombined_Sq = distUpperDepth_Sq + distAxis*distAxis;

#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2 )
        {
          G4double distCombined= std::sqrt( distCombined_Sq );
          G4double distUpperDepth= std::sqrt ( distUpperDepth_Sq );
          G4cout << " -- It is a Header " << G4endl;
          G4cout << "  Recurse to deal with next level, fVoxelDepth= " 
                 << fVoxelDepth+1 << G4endl;
          G4cout << "  Distances:  Upper= " << distUpperDepth
                 << " Axis= " << distAxis
                 << " Combined= " << distCombined << G4endl;
        }
#endif
        
        // Recurse to deal with lower levels
        //
        headerSafety= SafetyForVoxelHeader( pNewVoxelHeader, localPoint,
                                            maxLength, currentPhysical,
                                            distCombined_Sq, minSafety);
        ourSafety = std::min( ourSafety, headerSafety );
       
#ifdef G4DEBUG_NAVIGATION
        if( fVerbose > 2 )
        { 
          G4cout << "      >> Header safety = " << headerSafety
                 << " our level Saf = " << ourSafety << G4endl;
        }
#endif
        trialUp =   pNewVoxelHeader->GetMaxEquivalentSliceNo()+1;
        trialDown = pNewVoxelHeader->GetMinEquivalentSliceNo()-1;
     }
     minSafety = std::min( minSafety, ourSafety );

     // Find next closest Voxel
     //    - first try: by simple subtraction
     //    - later:  using distance  (TODO - tbc)
     //
     if( targetNodeNo >= pointNodeNo )
     {
        nextUp = trialUp;
        // distUp   = std::max( targetHeaderMax-localCrd,  0.0 );
        G4double lowerEdgeOfNext = targetHeaderMin
                                 + nextUp * targetHeaderNodeWidth;
        distUp = lowerEdgeOfNext-localCrd ;
        if( distUp < 0.0 )
        {
           distUp = DBL_MAX;  // On the wrong side - must not be considered
        }
#ifdef G4DEBUG_NAVIGATION
       if( fVerbose > 2 )
       {
         G4cout << " > Updated nextUp= " << nextUp << G4endl;
       }
#endif
     }
     
     if( targetNodeNo <= pointNodeNo ) 
     {
         nextDown = trialDown;
         // distDown = std::max( localCrd-targetHeaderMin,  0.0);
         G4double upperEdgeOfNext = targetHeaderMin
                                  + (nextDown+1) * targetHeaderNodeWidth;
         distDown = localCrd-upperEdgeOfNext;
         if( distDown < 0.0 )
         {
            distDown= DBL_MAX; // On the wrong side - must not be considered 
         }
#ifdef G4DEBUG_NAVIGATION
       if( fVerbose > 2 )
       {
         G4cout << " > Updated nextDown= " << nextDown << G4endl;
       }
#endif
     }

#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 2 )
     {
       G4cout << " Node= " << pointNodeNo
              << " Up:   next= " << nextUp  << " d# "
              << nextUp - pointNodeNo
              << " trialUp=  " << trialUp << " d# "
              << trialUp - pointNodeNo
              << " Down: next= " << nextDown << " d# "
              << targetNodeNo - nextDown
              << " trialDwn= " << trialDown << " d# "
              << targetNodeNo - trialDown
              << " condition (next is Inside)= " << nextIsInside
              << G4endl;
     }
#endif     
     
     G4bool UpIsClosest;
     UpIsClosest = distUp < distDown;
     
     if( (nextUp < targetHeaderNoSlices)
         && (UpIsClosest || (nextDown < 0)) )
     {
       nextNodeNo = nextUp;
       distAxis = distUp;
       ++nextUp; // Default
#ifdef G4VERBOSE
       if( fVerbose > 2 )
       {
          G4cout << " > Chose Up.   Depth= " << fVoxelDepth
                 << "   Nodes: next= " << nextNodeNo
                 << " new nextUp=   " << nextUp
                 << " Dist = " << distAxis << G4endl;
       }
#endif
     }
     else
     {
       nextNodeNo = nextDown;
       distAxis = distDown;
       --nextDown; // A default value
#ifdef G4VERBOSE
       if( fVerbose > 2 )
       {
         G4cout << " > Chose Down.  Depth= " << fVoxelDepth
                << "  Nodes: next= " << nextNodeNo
                << " new nextDown=  " << nextDown
                << " Dist = " << distAxis << G4endl;
       }
#endif
     }

     nextIsInside = (nextNodeNo >= 0) && (nextNodeNo < targetHeaderNoSlices);
     if( nextIsInside )
     {
       targetNodeNo= nextNodeNo;

#ifdef G4DEBUG_NAVIGATION
       assert( targetVoxelHeader->GetSlice(nextNodeNo) != 0 );
       G4bool bContinue = (distAxis<minSafety);
       if( !bContinue )
       {
         if( fVerbose > 2 )
         {
           G4cout << " Can skip remaining at depth " << targetHeaderAxis
                  << " >>  distAxis= " << distAxis
                  << " minSafety= " << minSafety << G4endl;
         }
       }
#endif
     }
     else
     {
#ifdef G4DEBUG_NAVIGATION
       if( fVerbose > 2)
       {
         G4cout << " VoxSaf> depth= " << fVoxelDepth << G4endl;
         G4cout << " VoxSaf> No more candidates:  nodeDown= " << nextDown
                << " nodeUp= " << nextUp
                << " lastSlice= " << targetHeaderNoSlices << G4endl;
       }
#endif
     }

     // This calculation can be 'hauled'-up to where minSafety is calculated
     //
     distMaxInterest = std::min( minSafety, maxLength ); 

  } while ( nextIsInside && ( distAxis*distAxis + distUpperDepth_Sq
                            < distMaxInterest*distMaxInterest ) ); 

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
  
  return ourSafety;
}
