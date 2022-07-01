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
// class G4SafetyHelper Implementation
//
// Original author: John Apostolakis, 2006
// --------------------------------------------------------------------

#include "G4SafetyHelper.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include "globals.hh"

G4SafetyHelper::G4SafetyHelper()
  : fLastSafetyPosition(0.0,0.0,0.0)
{
}

void G4SafetyHelper::InitialiseNavigator()
{
  fpPathFinder = G4PathFinder::GetInstance();
 
  G4TransportationManager* pTransportMgr= 
    G4TransportationManager::GetTransportationManager();

  fpMassNavigator = pTransportMgr->GetNavigatorForTracking(); 

  // Check
  //
  G4VPhysicalVolume* worldPV = fpMassNavigator->GetWorldVolume(); 
  if( worldPV == nullptr )
  { 
    G4Exception("G4SafetyHelper::InitialiseNavigator",
                "GeomNav0003", FatalException, 
                "Found that existing tracking Navigator has NULL world"); 
  }
}

void G4SafetyHelper::InitialiseHelper()
{
  fLastSafetyPosition = G4ThreeVector(0.0,0.0,0.0);
  fLastSafety         = 0.0;
  if (fFirstCall) { InitialiseNavigator(); }
  fFirstCall = false;
}

G4SafetyHelper::~G4SafetyHelper()
{
}

G4double   
G4SafetyHelper::CheckNextStep(const G4ThreeVector& position, 
                              const G4ThreeVector& direction,
                              const G4double currentMaxStep,
                                    G4double& newSafety )
{  
  // Distance in the Mass geometry
  //
  G4double linstep = fpMassNavigator->CheckNextStep( position,
                                                     direction,
                                                     currentMaxStep,
                                                     newSafety);
  fLastSafetyPosition = position;
  fLastSafety         = newSafety;

  // TO-DO: Can replace this with a call to PathFinder 
  //        giving id of Mass Geometry --> this avoid doing the work twice

  return linstep;
}

G4double G4SafetyHelper::ComputeSafety( const G4ThreeVector& position,
                                              G4double maxLength )
{
  G4double newSafety;
  
  // Only recompute (calling Navigator/PathFinder) if 'position'
  // is  *not* the safety location and has moved 'significantly'
  //
  G4double moveLengthSq = (position-fLastSafetyPosition).mag2();
  if(   (moveLengthSq > 0.0 ) )
  {
    if( !fUseParallelGeometries )
    {
      // Safety for mass geometry
      newSafety = fpMassNavigator->ComputeSafety(position, maxLength, true);

      // Only store a 'true' safety - one that was not restricted by maxLength
      if( newSafety < maxLength )
      {
         fLastSafety= newSafety;
         fLastSafetyPosition = position;
      }
    }
    else
    {
      // Safety for all geometries
      newSafety = fpPathFinder->ComputeSafety(position);

      fLastSafety= newSafety;
      fLastSafetyPosition = position;
    } 
 
  }
  else
  {
    // return last value if position is not (significantly) changed
    //
    // G4double moveLength = 0;
    // if( moveLengthSq > 0.0 ) { moveLength= std::sqrt(moveLengthSq); }
    newSafety = fLastSafety; // -moveLength;
  }
  
  return newSafety;
}

void G4SafetyHelper::ReLocateWithinVolume( const G4ThreeVector& newPosition )
{
#ifdef G4VERBOSE
  if( fVerbose > 0 )
  {
    // There is an opportunity - and need - to check whether
    // the proposed move is safe
    G4ThreeVector moveVec = newPosition - fLastSafetyPosition;
    if( moveVec.mag2() > sqr(fLastSafety) )
    {
      // A problem exists - we are proposing to move outside 'Safety Sphere'
      G4ExceptionDescription ed;
      ed << "Unsafe Move> Asked to relocate beyond 'Safety sphere'.  Details: "
         << G4endl;
      ed << " Safety Sphere:  Radius = " << fLastSafety;
      ed << " Center   = " << fLastSafetyPosition << G4endl;
      ed << " New Location :  Move   = " << moveVec.mag();
      ed << " Position = " << newPosition << G4endl;
      G4Exception("G4SafetyHelper::ReLocateWithinVolume",
                  "GeomNav1001", JustWarning, ed);
    }
  }
#endif
   
  if( !fUseParallelGeometries )
  {
    fpMassNavigator->LocateGlobalPointWithinVolume( newPosition ); 
  }
  else
  {
    fpPathFinder->ReLocate( newPosition ); 
  }
}

void  G4SafetyHelper::Locate( const G4ThreeVector& newPosition, 
                              const G4ThreeVector& newDirection)
{
  if( !fUseParallelGeometries )
  {
    fpMassNavigator->LocateGlobalPointAndSetup(newPosition, &newDirection,
                                               true, false); 
  }
  else
  {
    fpPathFinder->Locate( newPosition, newDirection ); 
  }
}
