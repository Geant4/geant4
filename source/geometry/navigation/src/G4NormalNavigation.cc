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
// $Id: G4NormalNavigation.cc,v 1.11 2010-11-04 08:57:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4NormalNavigation Implementation
//
// Author: P.Kent, 1996
//
// --------------------------------------------------------------------

#include "G4NormalNavigation.hh"
#include "G4AffineTransform.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4NormalNavigation::G4NormalNavigation()
  : fCheck(false)
{
  fLogger = new G4NavigationLogger("G4NormalNavigation");
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4NormalNavigation::~G4NormalNavigation()
{
  delete fLogger;
}

// ********************************************************************
// ComputeStep
// ********************************************************************
//
G4double
G4NormalNavigation::ComputeStep(const G4ThreeVector &localPoint,
                                const G4ThreeVector &localDirection,
                                const G4double currentProposedStepLength,
                                      G4double &newSafety,
                                      G4NavigationHistory &history,
                                      G4bool &validExitNormal,
                                      G4ThreeVector &exitNormal,
                                      G4bool &exiting,
                                      G4bool &entering,
                                      G4VPhysicalVolume *(*pBlockedPhysical),
                                      G4int &blockedReplicaNo)
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical, *blockedExitedVol=0;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep=currentProposedStepLength, motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;

  motherPhysical = history.GetTopVolume();
  motherLogical  = motherPhysical->GetLogicalVolume();
  motherSolid    = motherLogical->GetSolid();

  // Compute mother safety
  //
  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety; // Working isotropic safety
  
#ifdef G4VERBOSE
  if ( fCheck )
  {
    fLogger->PreComputeStepLog(motherPhysical, motherSafety, localPoint);
  }
#endif

  //
  // Compute daughter safeties & intersections
  //

  // Exiting normal optimisation
  //
  if ( exiting&&validExitNormal )
  {
    if ( localDirection.dot(exitNormal)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol =* pBlockedPhysical;
      ourSafety = 0;
    }
  }
  exiting  = false;
  entering = false;

  localNoDaughters = motherLogical->GetNoDaughters();
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo--)
  {
    samplePhysical = motherLogical->GetDaughter(sampleNo);
    if ( samplePhysical!=blockedExitedVol )
    {
      G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                 samplePhysical->GetTranslation());
      sampleTf.Invert();
      const G4ThreeVector samplePoint =
              sampleTf.TransformPoint(localPoint);
      const G4VSolid *sampleSolid =
              samplePhysical->GetLogicalVolume()->GetSolid();
      const G4double sampleSafety =
              sampleSolid->DistanceToIn(samplePoint);
#ifdef G4VERBOSE
      if( fCheck )
      {
        fLogger->PrintDaughterLog(sampleSolid, samplePoint, sampleSafety, 0);
      }
#endif
      if ( sampleSafety<ourSafety )
      {
        ourSafety=sampleSafety;
      }
      if ( sampleSafety<=ourStep )
      {
        sampleDirection = sampleTf.TransformAxis(localDirection);
        const G4double sampleStep =
                sampleSolid->DistanceToIn(samplePoint,sampleDirection);

#ifdef G4VERBOSE
        if( fCheck )
        {
          fLogger->PrintDaughterLog(sampleSolid, samplePoint,
                                    sampleSafety, sampleStep);
        }
#endif
        if ( sampleStep<=ourStep )
        {
          ourStep  = sampleStep;
          entering = true;
          exiting  = false;
          *pBlockedPhysical = samplePhysical;
          blockedReplicaNo  = -1;
#ifdef G4VERBOSE
          if( fCheck )
          {
            fLogger->AlongComputeStepLog(sampleSolid, samplePoint,
              sampleDirection, localDirection, sampleSafety, sampleStep);
          }
#endif
        }
      }
    }
  }
  if ( currentProposedStepLength<ourSafety )
  {
    // Guaranteed physics limited
    //
    entering = false;
    exiting  = false;
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
#ifdef G4VERBOSE
      if ( fCheck )
      {
        fLogger->PostComputeStepLog(motherSolid, localPoint, localDirection,
                                    motherStep, motherSafety);
      }
#endif

      if ( motherStep<=ourStep )
      {
        ourStep  = motherStep;
        exiting  = true;
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
  newSafety = ourSafety;
  return ourStep;
}

// ********************************************************************
// ComputeSafety
// ********************************************************************
//
G4double G4NormalNavigation::ComputeSafety(const G4ThreeVector &localPoint,
                                           const G4NavigationHistory &history,
                                           const G4double)
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4int localNoDaughters, sampleNo;

  motherPhysical = history.GetTopVolume();
  motherLogical  = motherPhysical->GetLogicalVolume();
  motherSolid    = motherLogical->GetSolid();

  // Compute mother safety
  //
  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety; // Working isotropic safety

#ifdef G4VERBOSE
  if( fCheck )
  {
    fLogger->ComputeSafetyLog(motherSolid, localPoint, motherSafety, true);
  }
#endif

  // Compute daughter safeties 
  //
  localNoDaughters = motherLogical->GetNoDaughters();
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical = motherLogical->GetDaughter(sampleNo);
    G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                               samplePhysical->GetTranslation());
    sampleTf.Invert();
    const G4ThreeVector samplePoint =
            sampleTf.TransformPoint(localPoint);
    const G4VSolid *sampleSolid =
            samplePhysical->GetLogicalVolume()->GetSolid();
    const G4double sampleSafety =
            sampleSolid->DistanceToIn(samplePoint);
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }
#ifdef G4VERBOSE
    if(fCheck)
    {
      fLogger->ComputeSafetyLog(sampleSolid, samplePoint, sampleSafety, false);
    }
#endif
  }
  return ourSafety;
}
