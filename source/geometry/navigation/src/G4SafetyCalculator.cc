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
// G4SafetyCalculator class Implementation
//
// Author: John Apostolakis, CERN - February 2023
// --------------------------------------------------------------------

#include "G4SafetyCalculator.hh"

#include "G4Navigator.hh"
#include "G4GeometryTolerance.hh"

G4SafetyCalculator::
G4SafetyCalculator( const G4Navigator& navigator,
                    const G4NavigationHistory& navHistory )
  : fNavigator(navigator), fNavHistory(navHistory)
{
  fkCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();   
}

// ********************************************************************
// ComputeSafety
//
// It assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the
//     ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.
// ********************************************************************
//
G4double G4SafetyCalculator::
SafetyInCurrentVolume( const G4ThreeVector& pGlobalpoint,
                             G4VPhysicalVolume* physicalVolume,
                       const G4double pMaxLength,
                             G4bool /* verbose */ )
{
  G4double safety = 0.0;
  G4ThreeVector stepEndPoint = fNavigator.GetLastStepEndPoint();

  G4ThreeVector localPoint = ComputeLocalPoint(pGlobalpoint);
  
  G4double distEndpointSq = (pGlobalpoint-stepEndPoint).mag2(); 
  G4bool stayedOnEndpoint = distEndpointSq < sqr(fkCarTolerance); 
  G4bool endpointOnSurface =   fNavigator.EnteredDaughterVolume()
                            || fNavigator.ExitedMotherVolume();
 
  G4VPhysicalVolume* motherPhysical = fNavHistory.GetTopVolume();
  if( motherPhysical != physicalVolume )
  {
     std::ostringstream msg;
     msg << " Current (navigation) phys-volume: " << motherPhysical
         << " name= " << motherPhysical->GetName() << G4endl
         << " Request made for     phys-volume: " << physicalVolume
         << " name= " << physicalVolume->GetName() << G4endl;
     G4Exception("G4SafetyCalculator::SafetyInCurrentVolume", "GeomNav0001",
                 FatalException, msg,
                 "This method must be called only in the Current volume.");
  }

  if( !(endpointOnSurface && stayedOnEndpoint) )
  {
    G4LogicalVolume* motherLogical = motherPhysical->GetLogicalVolume();
    G4SmartVoxelHeader* pVoxelHeader = motherLogical->GetVoxelHeader();

    // Pseudo-relocate to this point (updates voxel information only)
    //
    QuickLocateWithinVolume( localPoint, motherPhysical ); 
    //*********************

    // switch(CharacteriseDaughters(motherLogical))
    auto dtype= CharacteriseDaughters(motherLogical);
    switch(dtype)
    {
      case kNormal:
        if ( pVoxelHeader ) 
        {
          // New way: best safety
          safety = fVoxelSafety.ComputeSafety(localPoint,
                                               *motherPhysical, pMaxLength);
        }
        else
        {
          safety=fnormalNav.ComputeSafety(localPoint,fNavHistory,pMaxLength);
        }
        break;
      case kParameterised:
        if( GetDaughtersRegularStructureId(motherLogical) != 1 )
        {
          safety=fparamNav.ComputeSafety(localPoint,fNavHistory,pMaxLength);
        }
        else  // Regular structure
        {
          safety=fregularNav.ComputeSafety(localPoint,fNavHistory,pMaxLength);
        }
        break;
      case kReplica:
        safety = freplicaNav.ComputeSafety(pGlobalpoint, localPoint,
                                           fNavHistory, pMaxLength);
        break;
      case kExternal:
        safety = fpExternalNav->ComputeSafety(localPoint, fNavHistory,
                                              pMaxLength);
        break;
    }

    // Remember last safety origin & value
    //
    fPreviousSftOrigin = pGlobalpoint;
    fPreviousSafety = safety;
  }

  return safety;
}

// ********************************************************************
// QuickLocateWithinVolume
//
// -> the state information of this Navigator and its subNavigators
//    is updated in order to start the next step at pGlobalpoint
// -> no check is performed whether pGlobalpoint is inside the 
//    original volume (this must be the case).
//
// Note: a direction could be added to the arguments, to aid in future
//       optional checking (via the old code below, flagged by OLD_LOCATE). 
//       [ This would be done only in verbose mode ]
//
// Adapted simplied from G4Navigator::LocateGlobalPointWithinVolume()
// ********************************************************************
//
void G4SafetyCalculator::
QuickLocateWithinVolume( const G4ThreeVector& pointLocal,
                               G4VPhysicalVolume* motherPhysical )
{
  // For the case of Voxel (or Parameterised) volume the respective 
  // sub-navigator must be messaged to update its voxel information etc

  // Update the state of the Sub Navigators 
  // - in particular any voxel information they store/cache
  //
  G4LogicalVolume*    motherLogical  = motherPhysical->GetLogicalVolume();
  G4SmartVoxelHeader* pVoxelHeader   = motherLogical->GetVoxelHeader();

  switch( CharacteriseDaughters(motherLogical) )
  {
    case kNormal:
      if ( pVoxelHeader )
      {
        fvoxelNav.VoxelLocate( pVoxelHeader, pointLocal );
      }
      break;
    case kParameterised:
      if( GetDaughtersRegularStructureId(motherLogical) != 1 )
      {
        // Resets state & returns voxel node
        //
        fparamNav.ParamVoxelLocate( pVoxelHeader, pointLocal );
      }
      break;
    case kReplica:
      // Nothing to do
      break;
    case kExternal:
      fpExternalNav->RelocateWithinVolume( motherPhysical,
                                           pointLocal );
      break;
  }
}

// ********************************************************************
// Accessor for custom external navigation.
// ********************************************************************
G4VExternalNavigation* G4SafetyCalculator::GetExternalNavigation() const
{
  return fpExternalNav;   
}

// ********************************************************************
// Modifier for custom external navigation.
// ********************************************************************
void G4SafetyCalculator::SetExternalNavigation(G4VExternalNavigation* eNav)
{
  fpExternalNav = eNav;
}

// ********************************************************************
// CompareSafetyValues
// ********************************************************************
void G4SafetyCalculator::
CompareSafetyValues( G4double oldSafety,
                     G4double newValue,
                     G4VPhysicalVolume* motherPhysical,
               const G4ThreeVector& globalPoint,
                     G4bool keepState,
                     G4double maxLength,
                     G4bool enteredDaughterVol,
                     G4bool exitedMotherVol )
{
  constexpr G4double reportThreshold= 3.0e-14;
    // At least warn if rel-error exceeds it
  constexpr G4double  errorThreshold= 1.0e-08;
    // Fatal if relative error is larger
  constexpr G4double epsilonLen= 1.0e-20;
    // Baseline minimal value for divisor

  const G4double oldSafetyPlus = std::fabs(oldSafety)+epsilonLen;
  if( std::fabs( newValue - oldSafety) > reportThreshold * oldSafetyPlus )
  {
    G4ExceptionSeverity severity= FatalException;
    std::ostringstream msg;
    G4double diff= (newValue-oldSafety);
    G4double relativeDiff= diff / oldSafetyPlus;

    msg << " New (G4SafetyCalculator) value *disagrees* by relative diff " << relativeDiff
        << " in physical volume '" << motherPhysical->GetName() << "' "
        << "copy-no = " << motherPhysical->GetCopyNo();
    if( enteredDaughterVol ) { msg << "  ( Just Entered new daughter volume. ) "; }
    if( exitedMotherVol )    { msg << "  ( Just Exited previous volume. ) "; }
    msg << G4endl;
    msg << " Safeties:   old= "  << std::setprecision(12) << oldSafety
        << "   trial " << newValue
        << "  new-old= " << std::setprecision(7) << diff <<  G4endl;

    if( std::fabs(diff) < errorThreshold * ( std::fabs(oldSafety)+1.0e-20 ) )
    {
      msg << " (tiny difference) ";
      severity= JustWarning;
    }
    else
    {
      msg << " (real difference) ";
      severity= FatalException;

      // Extra information -- for big errors
      msg << " NOTE:  keepState =  " << keepState << G4endl;
      msg << " Location -  Global coordinates: " << globalPoint
          << "  volume= '" << motherPhysical->GetName() << "'"
          << " copy-no= " << motherPhysical->GetCopyNo() << G4endl;
      msg << " Argument maxLength= " << maxLength << G4endl;

      std::size_t depth= fNavHistory.GetDepth();
      msg << " Navigation History: depth = " << depth << G4endl;
      for( G4int i=1; i<(G4int)depth; ++i )
      {
         msg << "     d= " << i << " " << std::setw(32)
             << fNavHistory.GetVolume(i)->GetName()
             << "  copyNo= " << fNavHistory.GetReplicaNo(i);
         msg << G4endl;
      }
    }

#ifdef G4DEBUG_NAVIGATION
    G4double redo= SafetyInCurrentVolume(globalPoint, motherPhysical,
                                         maxLength, true);
    msg << " Redoing estimator: value = " << std::setprecision(16) << redo
        << "  diff/last= " << std::setprecision(7) << redo - newValue
        << "  diff/old= " << redo - oldSafety << G4endl;
#endif

    G4Exception("G4SafetyCalculator::CompareSafetyValues()", "GeomNav1007",
                severity, msg);
  }
}
