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
// class G4NormalNavigation Implementation
//
// Author: P.Kent, 1996
//
// --------------------------------------------------------------------

#include "G4NormalNavigation.hh"
#include "G4NavigationLogger.hh"
#include "G4AffineTransform.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4NormalNavigation::G4NormalNavigation()
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
//  On entry
//    exitNormal, validExitNormal:  for previous exited volume (daughter)
// 
//  On exit
//    exitNormal, validExitNormal:  for mother, if exiting it (else unchanged)
G4double
G4NormalNavigation::ComputeStep(const G4ThreeVector& localPoint,
                                const G4ThreeVector& localDirection,
                                const G4double currentProposedStepLength,
                                      G4double& newSafety,
                                      G4NavigationHistory& history,
                                      G4bool& validExitNormal,
                                      G4ThreeVector& exitNormal,
                                      G4bool& exiting,
                                      G4bool& entering,
                                      G4VPhysicalVolume* (*pBlockedPhysical),
                                      G4int& blockedReplicaNo)
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical,
                    *blockedExitedVol = nullptr;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4ThreeVector sampleDirection;
  G4double ourStep = currentProposedStepLength, ourSafety;
  G4double motherSafety, motherStep = DBL_MAX;
  G4long localNoDaughters, sampleNo;
  G4bool motherValidExitNormal = false;
  G4ThreeVector motherExitNormal; 

  motherPhysical = history.GetTopVolume();
  motherLogical  = motherPhysical->GetLogicalVolume();
  motherSolid    = motherLogical->GetSolid();

  // Compute mother safety
  //
  motherSafety = motherSolid->DistanceToOut(localPoint);
  ourSafety = motherSafety; // Working isotropic safety

  localNoDaughters = motherLogical->GetNoDaughters();
  
#ifdef G4VERBOSE
  if ( fCheck && ( (localNoDaughters>0) || (ourStep < motherSafety) )  )
  {
    fLogger->PreComputeStepLog(motherPhysical, motherSafety, localPoint);
  }
#endif
  // Compute daughter safeties & intersections
  //

  // Exiting normal optimisation
  //
  if ( exiting && validExitNormal )
  {
    if ( localDirection.dot(exitNormal)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol = (*pBlockedPhysical);
      ourSafety = 0;
    }
  }
  exiting  = false;
  entering = false;

#ifdef G4VERBOSE
  if ( fCheck )
  {
    // Compute early:
    //  a) to check whether point is (wrongly) outside
    //               (signaled if step < 0 or step == kInfinity )
    //  b) to check value against answer of daughters!

    motherStep = motherSolid->DistanceToOut(localPoint,
                                            localDirection,
                                            true,
                                           &motherValidExitNormal,
                                           &motherExitNormal);

    if( (motherStep >= kInfinity) || (motherStep < 0.0) )
    {
      // Error - indication of being outside solid !!
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

  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo--)
  {
    samplePhysical = motherLogical->GetDaughter(sampleNo);
    if ( samplePhysical!=blockedExitedVol )
    {
      G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                 samplePhysical->GetTranslation());
      sampleTf.Invert();
      const G4ThreeVector samplePoint = sampleTf.TransformPoint(localPoint);
      const G4VSolid *sampleSolid =
              samplePhysical->GetLogicalVolume()->GetSolid();
      const G4double sampleSafety =
              sampleSolid->DistanceToIn(samplePoint);

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
                                    sampleSafety, true,
                                    sampleDirection, sampleStep);          
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
                                         sampleDirection, localDirection,
                                         sampleSafety, sampleStep);
          }
#endif          
        }

#ifdef G4VERBOSE
        if( fCheck && (sampleStep < kInfinity) && (sampleStep >= motherStep) )
        {
           // The intersection point with the daughter is at or after the exit
           // point from the mother volume.  Double check!
           fLogger->CheckDaughterEntryPoint(sampleSolid,
                                            samplePoint, sampleDirection,
                                            motherSolid,
                                            localPoint, localDirection,
                                            motherStep, sampleStep);
        }
#endif
      } // end of if ( sampleSafety <= ourStep ) 
#ifdef G4VERBOSE
      else if ( fCheck )
      {
         fLogger->PrintDaughterLog(sampleSolid, samplePoint,
                                   sampleSafety, false,
                                   G4ThreeVector(0.,0.,0.), -1.0 );
      }
#endif          
    }
  }
  if ( currentProposedStepLength<ourSafety )
  {
    // Guaranteed physics limited
    //
    entering = false;
    exiting  = false;
    *pBlockedPhysical = nullptr;
    ourStep = kInfinity;
  }
  else
  {
    // Consider intersection with mother solid
    //
    if ( motherSafety<=ourStep )
    {
      if ( !fCheck )  // The call is moved above when running in check_mode
      {
        motherStep = motherSolid->DistanceToOut(localPoint,
                                                localDirection,
                                                true,
                                               &motherValidExitNormal,
                                               &motherExitNormal);
      }
#ifdef G4VERBOSE
      else  // check_mode
      {
        fLogger->PostComputeStepLog(motherSolid, localPoint, localDirection,
                                    motherStep, motherSafety);
        if( motherValidExitNormal )
        {
          fLogger->CheckAndReportBadNormal(motherExitNormal,
                                           localPoint,
                                           localDirection,
                                           motherStep,
                                           motherSolid,
                                           "From motherSolid::DistanceToOut" );
        }
      }
#endif

      if( (motherStep >= kInfinity) || (motherStep < 0.0) )
      {
#ifdef G4VERBOSE
        if( fCheck )  // Clearly outside the mother solid!
        {
          fLogger->ReportOutsideMother(localPoint, localDirection,
                                       motherPhysical);
        }
#endif
        ourStep = motherStep = 0.0;
        exiting = true;
        entering = false;
        // validExitNormal= motherValidExitNormal;
        // exitNormal= motherExitNormal;
        //  The normal could be useful - but only if near the mother
        //  But it could be unreliable!
        validExitNormal = false;
        *pBlockedPhysical = nullptr; // or motherPhysical ?
        blockedReplicaNo = 0;  // or motherReplicaNumber ?
        newSafety= 0.0;
        return ourStep;
      }

      if ( motherStep<=ourStep )
      {
        ourStep  = motherStep;
        exiting  = true;
        entering = false;
        validExitNormal = motherValidExitNormal;
        exitNormal = motherExitNormal;
        
        if ( motherValidExitNormal )
        {
          const G4RotationMatrix *rot = motherPhysical->GetRotation();
          if (rot)
          {
            exitNormal *= rot->inverse();
#ifdef G4VERBOSE
            if( fCheck )
               fLogger->CheckAndReportBadNormal(exitNormal,        // rotated
                                                motherExitNormal,  // original 
                                                *rot,
                                                "From RotationMatrix" );
#endif            
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
G4double G4NormalNavigation::ComputeSafety(const G4ThreeVector& localPoint,
                                           const G4NavigationHistory& history,
                                           const G4double)
{
  G4VPhysicalVolume *motherPhysical, *samplePhysical;
  G4LogicalVolume *motherLogical;
  G4VSolid *motherSolid;
  G4double motherSafety, ourSafety;
  G4long localNoDaughters, sampleNo;

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
    fLogger->ComputeSafetyLog(motherSolid,localPoint,motherSafety,true,true);
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
      fLogger->ComputeSafetyLog(sampleSolid, samplePoint,
                                sampleSafety, false, false);
        // Not mother, no banner
    }
#endif
  }
  return ourSafety;
}

// The following methods have been imported to this source file
// in order to avoid dependency of the header file on the
// header implementation of G4NavigationLogger.

// ********************************************************************
// GetVerboseLevel
// ********************************************************************
//
G4int G4NormalNavigation::GetVerboseLevel() const
{
  return fLogger->GetVerboseLevel();
}

// ********************************************************************
// SetVerboseLevel
// ********************************************************************
//
void G4NormalNavigation::SetVerboseLevel(G4int level)
{
  fLogger->SetVerboseLevel(level);
}
