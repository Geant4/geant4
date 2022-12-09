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
// class G4MultiNavigator Implementation
//
// Author:  John Apostolakis, November 2006
// --------------------------------------------------------------------

#include <iomanip>

#include "G4MultiNavigator.hh"

class G4FieldManager;

#include "G4SystemOfUnits.hh"
#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4MultiNavigator::G4MultiNavigator() 
  : G4Navigator()
{
  G4ThreeVector Big3Vector( kInfinity, kInfinity, kInfinity ); 
  fLastLocatedPosition = Big3Vector;
  fSafetyLocation  = Big3Vector;
  fPreStepLocation = Big3Vector;

  for(auto num=0; num< fMaxNav; ++num )
  {
    fpNavigator[num] = nullptr;   
    fLimitTruth[num] = false;
    fLimitedStep[num] = kUndefLimited;
    fCurrentStepSize[num] = fNewSafety[num] = -1.0; 
    fLocatedVolume[num] = nullptr; 
  }

  pTransportManager= G4TransportationManager::GetTransportationManager();

  G4Navigator* massNav= pTransportManager->GetNavigatorForTracking();
  if( massNav )
  { 
    G4VPhysicalVolume* pWorld= massNav->GetWorldVolume(); 
    if( pWorld )
    { 
      SetWorldVolume( pWorld ); 
      fLastMassWorld = pWorld; 
    }
  }
}

G4MultiNavigator::~G4MultiNavigator() 
{
}

G4double G4MultiNavigator::ComputeStep(const G4ThreeVector& pGlobalPoint,
                                       const G4ThreeVector& pDirection,
                                       const G4double       proposedStepLength,
                                             G4double&      pNewSafety)
{
  G4double safety= 0.0, step=0.0;
  G4double minSafety= kInfinity, minStep= kInfinity;

  fNoLimitingStep= -1; 
  fIdNavLimiting= -1;     // Reset for new step

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 2 )
  {
    G4cout << " G4MultiNavigator::ComputeStep : entered " << G4endl;
    G4cout << "   Input position= " << pGlobalPoint
           << "   direction= "      << pDirection         << G4endl;
    G4cout << "   Requested step= " << proposedStepLength << G4endl;
  }
#endif

  std::vector<G4Navigator*>::iterator pNavigatorIter;

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

  G4ThreeVector initialPosition = pGlobalPoint;
  G4ThreeVector initialDirection= pDirection;

  for( auto num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num )
  {
     safety= kInfinity;

     step= (*pNavigatorIter)->ComputeStep( initialPosition, 
                                           initialDirection,
                                           proposedStepLength,
                                           safety ); 
     if( safety < minSafety ){ minSafety = safety; } 
     if( step < minStep )    { minStep= step; } 

     fCurrentStepSize[num] = step; 
     fNewSafety[num]= safety; 
      // This is currently the safety from the last sub-step

#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 2 )
     {
       G4cout << "G4MultiNavigator::ComputeStep : Navigator ["
              << num << "] -- step size " << step
              << " safety= " << safety << G4endl;
     }
#endif
  } 

  // Save safety value, related position
  //
  fPreStepLocation     = initialPosition; 
  fMinSafety_PreStepPt = minSafety;
  fMinStep = minStep; 

  if( fMinStep == kInfinity )
  {
     fTrueMinStep = proposedStepLength;   //  Use this below for endpoint !!
  }
  else
  {
     fTrueMinStep = minStep;
  }

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 1 )
  {
    G4ThreeVector endPosition = initialPosition+fTrueMinStep*initialDirection;

    G4int oldPrec = G4cout.precision(8); 
    G4cout << "G4MultiNavigator::ComputeStep : "
           << " initialPosition = " << initialPosition 
           << " and endPosition = " << endPosition << G4endl;
    G4cout.precision( oldPrec );
  }
#endif

  pNewSafety = minSafety; 

  this->WhichLimited(); 

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 2 )
  {
    G4cout << " G4MultiNavigator::ComputeStep : exits returning "
           << minStep << G4endl;
  }
#endif

  return minStep;  // must return kInfinity if do not limit step
}

// ----------------------------------------------------------------------

G4double 
G4MultiNavigator::ObtainFinalStep( G4int     navigatorId, 
                                   G4double& pNewSafety,  // for this geometry
                                   G4double& minStep,
                                   ELimited& limitedStep) 
{
  if( navigatorId > fNoActiveNavigators )
  { 
     std::ostringstream message;
     message << "Bad Navigator Id!" << G4endl
             << "        Navigator Id = " << navigatorId 
             << "        No Active = " << fNoActiveNavigators << ".";
     G4Exception("G4MultiNavigator::ObtainFinalStep()", "GeomNav0002",
                 FatalException, message); 
  }

  // Prepare the information to return
  //
  pNewSafety  = fNewSafety[ navigatorId ]; 
  limitedStep = fLimitedStep[ navigatorId ];
  minStep= fMinStep; 

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 1 )
  { 
     G4cout << " G4MultiNavigator::ComputeStep returns "
            << fCurrentStepSize[ navigatorId ]
            << " for Navigator " << navigatorId
            << " Limited step = " << limitedStep 
            << " Safety(mm) = " << pNewSafety / mm << G4endl; 
  }
#endif

  return fCurrentStepSize[ navigatorId ];
}

// ----------------------------------------------------------------------

void G4MultiNavigator::PrepareNewTrack( const G4ThreeVector position, 
                                        const G4ThreeVector direction )
{
#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 1 )
  {
    G4cout << " Entered G4MultiNavigator::PrepareNewTrack() " << G4endl;
  }
#endif

  G4MultiNavigator::PrepareNavigators(); 

  LocateGlobalPointAndSetup( position, &direction, false, false );   
  //
  // The first location for each Navigator must be non-relative
  // or else call ResetStackAndState() for each Navigator
  // Use direction to get correct side of boundary (ignore dir= false)
}

// ----------------------------------------------------------------------

void G4MultiNavigator::PrepareNavigators()
{
  // Key purposes:
  //   - Check and cache set of active navigators
  //   - Reset state for new track

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 1 )
  {
    G4cout << " Entered G4MultiNavigator::PrepareNavigators() " << G4endl;
  }
#endif

  // Message the transportation-manager to find active navigators

  std::vector<G4Navigator*>::const_iterator pNavigatorIter; 
  fNoActiveNavigators = (G4int)pTransportManager-> GetNoActiveNavigators();

  if( fNoActiveNavigators > fMaxNav )
  {
    std::ostringstream message;
    message << "Too many active Navigators / worlds !" << G4endl
            << "        Active Navigators (worlds): "
            << fNoActiveNavigators << G4endl
            << "        which is more than the number allowed: "
            << fMaxNav << " !";
    G4Exception("G4MultiNavigator::PrepareNavigators()", "GeomNav0002",  
                FatalException, message); 
  }

  pNavigatorIter= pTransportManager->GetActiveNavigatorsIterator();
  for( auto num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num )
  {
     fpNavigator[num] = *pNavigatorIter;   
     fLimitTruth[num] = false;
     fLimitedStep[num] = kDoNot;
     fCurrentStepSize[num] = 0.0; 
     fLocatedVolume[num] = nullptr; 
  }
  fWasLimitedByGeometry = false; 

  // Check the world volume of the mass navigator
  // in case a call to SetWorldVolume() changed it

  G4VPhysicalVolume* massWorld = GetWorldVolume();

  if( (massWorld != fLastMassWorld) && (massWorld!=0) )
  { 
     // Pass along change to Mass Navigator
     fpNavigator[0] -> SetWorldVolume( massWorld );

#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 0 )
     { 
       G4cout << " G4MultiNavigator::PrepareNavigators() changed world volume " 
              << " for mass geometry to " << massWorld->GetName() << G4endl; 
     }
#endif

     fLastMassWorld = massWorld;
  }
}

// ----------------------------------------------------------------------

G4VPhysicalVolume* 
G4MultiNavigator::LocateGlobalPointAndSetup(const G4ThreeVector& position,
                                            const G4ThreeVector* pDirection,
                                            const G4bool pRelativeSearch,
                                            const G4bool ignoreDirection )
{
  // Locate the point in each geometry

  G4ThreeVector direction(0.0, 0.0, 0.0);
  G4bool relative = pRelativeSearch; 
  std::vector<G4Navigator*>::iterator pNavIter
    = pTransportManager->GetActiveNavigatorsIterator(); 

  if( pDirection ) { direction = *pDirection; }

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 2 )
  {
    G4cout << " Entered G4MultiNavigator::LocateGlobalPointAndSetup() "
           << G4endl;
    G4cout << "   Locating at position: " << position
           << ", with direction: " << direction << G4endl
           << "   Relative: " << relative 
           << ", ignore direction: " << ignoreDirection << G4endl;
    G4cout << "   Number of active navigators: " << fNoActiveNavigators
           << G4endl;
  }
#endif

  for ( auto num=0; num< fNoActiveNavigators ; ++pNavIter,++num )
  {
     if( fWasLimitedByGeometry && fLimitTruth[num] )
     { 
        (*pNavIter)->SetGeometricallyLimitedStep(); 
     }

     G4VPhysicalVolume *pLocated
       = (*pNavIter)->LocateGlobalPointAndSetup( position, &direction,
                                                 relative, ignoreDirection );   
     // Set the state related to the location
     //
     fLocatedVolume[num] = pLocated; 

     // Clear state related to the step
     //
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
     fLimitTruth[ num ] = false;   // Always clear on locating (see Navigator)
    
#ifdef G4DEBUG_NAVIGATION
     if( fVerbose > 2 )
     {
       G4cout << " Located in world: " << num << ", at: " << position << G4endl
              << " Used geomLimStp: " << fLimitTruth[num]
              << ", found in volume: " << pLocated << G4endl; 
       G4cout << " Name = '" ;       
       if( pLocated )
       { 
         G4cout << pLocated->GetName() << "'"; 
         G4cout << " - CopyNo= " << pLocated->GetCopyNo(); 
       }
       else
       { 
         G4cout <<  "Null'   Id: Not-Set "; 
       }
       G4cout << G4endl; 
     }
#endif
  }

  fWasLimitedByGeometry = false;   // Clear on locating
  G4VPhysicalVolume* volMassLocated = fLocatedVolume[0]; 

  return volMassLocated;
}

// ----------------------------------------------------------------------

void
G4MultiNavigator::LocateGlobalPointWithinVolume(const G4ThreeVector& position)
{
  // Relocate the point in each geometry

  std::vector<G4Navigator*>::iterator pNavIter
    = pTransportManager->GetActiveNavigatorsIterator(); 

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 2 )
  {
    G4cout << " Entered G4MultiNavigator::ReLocate() " << G4endl
           << "  Re-locating at position: " << position  << G4endl; 
  }
#endif

  for ( auto num=0; num< fNoActiveNavigators ; ++pNavIter,++num )
  {
     //  ... none limited the step

     (*pNavIter)->LocateGlobalPointWithinVolume( position ); 

     // Clear state related to the step
     //
     fLimitedStep[num]     = kDoNot; 
     fCurrentStepSize[num] = 0.0;      

     fLimitTruth[ num ] = false;   // Always clear on locating (see Navigator)
  }
  fWasLimitedByGeometry = false;   // Clear on locating
  fLastLocatedPosition  = position; 
}

// ----------------------------------------------------------------------

G4double G4MultiNavigator::ComputeSafety( const G4ThreeVector& position,
                                          const G4double       maxDistance,
                                          const G4bool         state)
{
    // Recompute safety for the relevant point

    G4double minSafety = kInfinity, safety = kInfinity;
  
    std::vector<G4Navigator*>::iterator pNavigatorIter;
    pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

    for( auto num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num )
    {
       safety = (*pNavigatorIter)->ComputeSafety( position, maxDistance, state);
       if( safety < minSafety ) { minSafety = safety; } 
    } 

    fSafetyLocation = position;
    fMinSafety_atSafLocation = minSafety;

#ifdef G4DEBUG_NAVIGATION
    if( fVerbose > 1 )
    { 
      G4cout << " G4MultiNavigator::ComputeSafety - returns: " 
             << minSafety << ", at location: " << position << G4endl;
    }
#endif
    return minSafety; 
}

// -----------------------------------------------------------------------

G4TouchableHistoryHandle 
G4MultiNavigator::CreateTouchableHistoryHandle() const
{
  G4Exception( "G4MultiNavigator::CreateTouchableHistoryHandle()", 
               "GeomNav0001", FatalException,  
               "Getting a touchable from G4MultiNavigator is not defined."); 

  G4TouchableHistory* touchHist;
  touchHist= fpNavigator[0] -> CreateTouchableHistory(); 

  G4VPhysicalVolume* locatedVolume= fLocatedVolume[0]; 
  if( locatedVolume == nullptr )
  {
    // Workaround to ensure that the touchable is fixed !! // TODO: fix
    //
    touchHist->UpdateYourself( locatedVolume, touchHist->GetHistory() );
  }
    
  return G4TouchableHistoryHandle(touchHist); 
}

// -----------------------------------------------------------------------

void G4MultiNavigator::WhichLimited()
{
  // Flag which processes limited the step

  G4int last=-1; 
  const G4int IdTransport= 0;  // Id of Mass Navigator !!
  G4int noLimited=0; 
  ELimited shared= kSharedOther; 

#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 2 )
  {
    G4cout << " Entered G4MultiNavigator::WhichLimited() " << G4endl;
  }
#endif

  // Assume that [IdTransport] is Mass / Transport
  // 
  G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep)
                         && ( fMinStep!= kInfinity); 
  if( transportLimited )
  { 
     shared = kSharedTransport;
  }

  for ( auto num = 0; num < fNoActiveNavigators; ++num )
  { 
    G4bool limitedStep;

    G4double step = fCurrentStepSize[num]; 

    limitedStep = ( step == fMinStep ) && ( step != kInfinity); 
   
    fLimitTruth[ num ] = limitedStep; 
    if( limitedStep )
    {
      ++noLimited;  
      fLimitedStep[num] = shared;
      last= num; 
    }
    else
    {
      fLimitedStep[num] = kDoNot;
    }
  }
  if( (last > -1) && (noLimited == 1 ) )
  {
    fLimitedStep[ last ] = kUnique;
    fIdNavLimiting = last;
  }

  fNoLimitingStep = noLimited;  

  return;
}

// -----------------------------------------------------------------------

void
G4MultiNavigator::PrintLimited()
{
  // Report results -- for checking   

  static const G4String StrDoNot("DoNot"), StrUnique("Unique"),
                  StrUndefined("Undefined"),
                  StrSharedTransport("SharedTransport"),
                  StrSharedOther("SharedOther");
  G4cout << "### G4MultiNavigator::PrintLimited() reports: " << G4endl;
  G4cout << "    Minimum step (true): " << fTrueMinStep 
         << ", reported min: " << fMinStep << G4endl; 

#ifdef G4DEBUG_NAVIGATION
  if(fVerbose>=2)
  {
    G4cout << std::setw(5) << " NavId"  << " "
           << std::setw(12) << " step-size " << " "
           << std::setw(12) << " raw-size "  << " "
           << std::setw(12) << " pre-safety " << " " 
           << std::setw(15) << " Limited / flag"  << " "
           << std::setw(15) << "  World "  << " "
           << G4endl;  
  }
#endif

  for ( auto num = 0; num < fNoActiveNavigators; ++num )
  { 
    G4double rawStep = fCurrentStepSize[num]; 
    G4double stepLen = fCurrentStepSize[num]; 
    if( stepLen > fTrueMinStep )
    { 
      stepLen = fTrueMinStep;     // did not limit (went as far as asked)
    }
    G4long oldPrec = G4cout.precision(9); 

    G4cout << std::setw(5) << num  << " "
           << std::setw(12) << stepLen << " "
           << std::setw(12) << rawStep << " "
           << std::setw(12) << fNewSafety[num] << " "
           << std::setw(5) << (fLimitTruth[num] ? "YES" : " NO") << " ";
    G4String limitedStr;
    switch ( fLimitedStep[num] )
    {
      case kDoNot          : limitedStr = StrDoNot; break;
      case kUnique         : limitedStr = StrUnique; break; 
      case kSharedTransport: limitedStr = StrSharedTransport; break; 
      case kSharedOther    : limitedStr = StrSharedOther; break;
      default              : limitedStr = StrUndefined; break;
    }
    G4cout << " " << std::setw(15) << limitedStr << " ";  
    G4cout.precision(oldPrec); 

    G4Navigator *pNav = fpNavigator[ num ]; 
    G4String  WorldName( "Not-Set" ); 
    if (pNav != nullptr)
    {
       G4VPhysicalVolume *pWorld = pNav->GetWorldVolume(); 
       if( pWorld )
       {
           WorldName = pWorld->GetName(); 
       }
    }
    G4cout << " " << WorldName ; 
    G4cout << G4endl;
  }
}
 

// -----------------------------------------------------------------------

void G4MultiNavigator::ResetState()
{
   fWasLimitedByGeometry = false; 

   G4Exception("G4MultiNavigator::ResetState()", "GeomNav0001",  
               FatalException,  
               "Cannot reset state for navigators of G4MultiNavigator.");
/*
   std::vector<G4Navigator*>::iterator pNavigatorIter;
   pNavigatorIter = pTransportManager->GetActiveNavigatorsIterator();
   for( auto num = 0; num< fNoActiveNavigators; ++pNavigatorIter,++num )
   {
       //  (*pNavigatorIter)->ResetState();  // KEEP THIS comment !!!
   }
*/
}

// -----------------------------------------------------------------------

void G4MultiNavigator::SetupHierarchy()
{
  G4Exception( "G4MultiNavigator::SetupHierarchy()", 
               "GeomNav0001", FatalException,  
               "Cannot setup hierarchy for navigators of G4MultiNavigator."); 
}

// -----------------------------------------------------------------------

void G4MultiNavigator::CheckMassWorld()
{
   G4VPhysicalVolume* navTrackWorld =
     pTransportManager->GetNavigatorForTracking()->GetWorldVolume();

   if( navTrackWorld != fLastMassWorld )
   { 
      G4Exception( "G4MultiNavigator::CheckMassWorld()",
                   "GeomNav0003", FatalException, 
                   "Mass world pointer has been changed." ); 
   }
}

// -----------------------------------------------------------------------

G4VPhysicalVolume*
G4MultiNavigator::ResetHierarchyAndLocate(const G4ThreeVector& point,
                                          const G4ThreeVector& direction,
                                          const G4TouchableHistory& MassHistory)
{
   // Reset geometry for all -- and use the touchable for the mass history

   G4VPhysicalVolume* massVolume = nullptr; 
   G4Navigator* pMassNavigator = fpNavigator[0]; 

   if( pMassNavigator != nullptr )
   {
      massVolume= pMassNavigator->ResetHierarchyAndLocate( point, direction,
                                                           MassHistory); 
   }
   else
   {
      G4Exception("G4MultiNavigator::ResetHierarchyAndLocate()",
                  "GeomNav0002", FatalException,
                  "Cannot reset hierarchy before navigators are initialised.");
   }

   std::vector<G4Navigator*>::iterator pNavIter= 
       pTransportManager->GetActiveNavigatorsIterator(); 

   for ( auto num = 0; num < fNoActiveNavigators ; ++pNavIter,++num )
   {
      G4bool relativeSearch, ignoreDirection; 

      (*pNavIter)-> LocateGlobalPointAndSetup( point, 
                                               &direction, 
                                               relativeSearch=false,
                                               ignoreDirection=false);
   }
   return massVolume; 
}

// -----------------------------------------------------------------------

G4ThreeVector 
G4MultiNavigator::GetGlobalExitNormal(const G4ThreeVector& argPoint,
                                      G4bool* argpObtained)  //  obtained valid
{
  G4ThreeVector normalGlobalCrd(0.0, 0.0, 0.0); 
  G4bool isObtained = false; 
  // These default values will be used if fNoLimitingStep== 0
  G4int firstNavigatorId = -1;
  G4bool oneObtained = false;
    
  if( fNoLimitingStep == 1 )
  { 
    // Only message the Navigator which limited the step!
    normalGlobalCrd = fpNavigator[ fIdNavLimiting ]
                    ->GetGlobalExitNormal( argPoint, &isObtained ); 
    *argpObtained = isObtained; 
  }
  else 
  {
    if( fNoLimitingStep > 1 )
    { 
      std::vector<G4Navigator*>::iterator pNavIter= 
        pTransportManager->GetActiveNavigatorsIterator(); 

      for ( auto num = 0; num < fNoActiveNavigators ; ++pNavIter, ++num )
      {
        G4ThreeVector oneNormal;
        if( fLimitTruth[ num ] )  // Did this geometry limit the step ? 
        { 
          G4ThreeVector newNormal =
            (*pNavIter)->GetGlobalExitNormal( argPoint, &oneObtained );
          if( oneObtained )
          {
            // Keep first one - only if it is valid (ie not null)
            if( !isObtained && (newNormal.mag2() != 0.0) )
            {
              normalGlobalCrd = newNormal; 
              isObtained = oneObtained;
              firstNavigatorId = num; 
            }
            else
            {
              // Check for clash
              G4double dotNewPrevious = newNormal.dot( normalGlobalCrd );
              G4double productMagSq = normalGlobalCrd.mag2()*newNormal.mag2(); 
              if( productMagSq > 0.0 ) 
              {
                G4double productMag = std::sqrt( productMagSq ); 
                dotNewPrevious /= productMag; // Normalise
                if( dotNewPrevious < (1 - perThousand) ) 
                {
                  *argpObtained = false; 

                  if( fVerbose > 2 )  //  dotNewPrevious <= 0.0 )
                  {
                    std::ostringstream message;
                    message << "Clash of Normal from different Navigators!"
                            << G4endl
                            << "        Previous Navigator Id = "
                            << firstNavigatorId << G4endl
                            << "        Current  Navigator Id = "
                            << num << G4endl;
                    message << "  Dot product of 2 normals = "
                            << dotNewPrevious << G4endl;
                    message << "        Normal (previous) = "
                            << normalGlobalCrd << G4endl;
                    message << "        Normal (current)  = " << newNormal       << G4endl;
                    G4Exception("G4MultiNavigator::GetGlobalExitNormal()",
                                "GeomNav0002", JustWarning, message); 
                  }
                }
                else
                {
                  // Close agreement - Do not change 
                }
              }
            }
          }
        }
      } // end for over the Navigators

      // Report if no Normal was obtained
      if( !oneObtained ) 
      {
        std::ostringstream message;
        message << "No Normal obtained despite having " << fNoLimitingStep
                << " candidate Navigators limiting the step!" << G4endl;
        G4Exception("G4MultiNavigator::GetGlobalExitNormal()", "GeomNav0002",
                    JustWarning, message); 
      }

    } // end if ( fNoLimiting > 1 ) 
  } // end else

  *argpObtained = isObtained;
  return normalGlobalCrd;
}

// -----------------------------------------------------------------------

G4ThreeVector 
G4MultiNavigator::GetLocalExitNormal(G4bool* argpObtained)
{
  // If it is the mass navigator, then expect
  G4ThreeVector normalGlobalCrd(0.0, 0.0, 0.0); 
  G4bool isObtained = false; 
  // These default values will be used if fNoLimitingStep== 0

  if( fNoLimitingStep == 1 )
  { 
    // Only message the Navigator which limited the step!
    normalGlobalCrd = fpNavigator[ fIdNavLimiting ]
                    ->GetLocalExitNormal( &isObtained ); 
    *argpObtained = isObtained;

    static G4ThreadLocal G4int numberWarnings = 0;
    G4int noWarningsStart = 10, noModuloWarnings = 100; 
    ++numberWarnings; 
    if( (numberWarnings < noWarningsStart )
     || (numberWarnings%noModuloWarnings == 0) ) 
    {
      std::ostringstream message;
      message << "Cannot obtain normal in local coordinates of two or more "
              << "coordinate systems." << G4endl;
      G4Exception("G4MultiNavigator::GetGlobalExitNormal()", "GeomNav0002",
                  JustWarning, message);       
    }
  }
  else
  {
    if( fNoLimitingStep > 1 ) 
    {
      std::ostringstream message;
      message << "Cannot obtain normal in local coordinates of two or more "
              << "coordinate systems." << G4endl;
      G4Exception("G4MultiNavigator::GetGlobalExitNormal()", "GeomNav0002",
                  FatalException, message);       
    }
  }
    
  *argpObtained = isObtained;
  return normalGlobalCrd; 
}

// -----------------------------------------------------------------------

G4ThreeVector 
G4MultiNavigator::GetLocalExitNormalAndCheck(const G4ThreeVector&, // point,
                                                   G4bool* obtained)
{
  return G4MultiNavigator::GetLocalExitNormal( obtained );
}
