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
// $Id: G4NormalNavigation.cc,v 1.10 2010-11-02 13:41:29 gcosmo Exp $
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

#include <iomanip>

// ********************************************************************
// Constructor
// ********************************************************************
//
G4NormalNavigation::G4NormalNavigation()
  : fCheck(false), fVerbose(0)
{
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4NormalNavigation::~G4NormalNavigation()
{;}

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
    PreComputeStepLog(motherPhysical, motherSafety, localPoint);
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
      if ( sampleSafety<ourSafety )
      {
        ourSafety=sampleSafety;
      }
      if ( sampleSafety<=ourStep )
      {
        sampleDirection = sampleTf.TransformAxis(localDirection);
        const G4double sampleStep =
                sampleSolid->DistanceToIn(samplePoint,sampleDirection);

        if ( sampleStep<=ourStep )
        {
          ourStep  = sampleStep;
          entering = true;
          exiting  = false;
          *pBlockedPhysical = samplePhysical;
          blockedReplicaNo  = -1;
        }
#ifdef G4VERBOSE
        if(( fCheck ) && ( fVerbose == 1 ))
        {
          AlongComputeStepLog(sampleSolid, samplePoint, sampleDirection,
                           localDirection, sampleSafety, sampleStep, ourStep);
        }
#endif
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
      PostComputeStepLog(motherSolid, localPoint, localDirection,
                         motherStep, motherSafety);
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
  if(( fCheck ) && ( fVerbose == 1 ))
  {
    ComputeSafetyLog(motherSolid, localPoint, motherSafety, true);
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
    if(( fCheck ) && ( fVerbose == 1 ))
    {
      ComputeSafetyLog(sampleSolid, samplePoint, sampleSafety, false);
    }
#endif
  }
  return ourSafety;
}

// ********************************************************************
// Logging methods ...
// ********************************************************************

void
G4NormalNavigation::PreComputeStepLog(const G4VPhysicalVolume* motherPhysical,
                                            G4double motherSafety,
                                      const G4ThreeVector& localPoint)
{
    G4VSolid* motherSolid = motherPhysical->GetLogicalVolume()->GetSolid();
    if( fVerbose == 1 )
    {
      G4cout << "*** G4NormalNavigation::ComputeStep(): ***" << G4endl
             << "    Invoked DistanceToOut(p) for mother solid: "
             << motherSolid->GetName()
             << ". Solid replied: " << motherSafety << G4endl
             << "    For local point p: " << localPoint << G4endl
             << "    To be considered as 'mother safety'." << G4endl;
    }
    if ( motherSafety < 0.0 )
    {
      G4cerr << "ERROR - G4NormalNavigation::ComputeStep()" << G4endl
             << "        Current solid " << motherSolid->GetName()
             << " gave negative safety: " << motherSafety << G4endl
             << "        for the current (local) point " << localPoint
             << G4endl;
      motherSolid->DumpInfo();
      G4Exception("G4NormalNavigation::ComputeStep()",
                  "NegativeSafetyMotherVol", FatalException,
                  "Negative Safety In Voxel Navigation !" );
    }
    if( motherSolid->Inside(localPoint)==kOutside )
    {
      G4cout << "WARNING - G4NormalNavigation::ComputeStep()" << G4endl
             << "          Point " << localPoint
             << " is outside current volume " << motherPhysical->GetName()
             << G4endl;
      G4double estDistToSolid= motherSolid->DistanceToIn(localPoint); 
      G4cout << "          Estimated isotropic distance to solid (distToIn)= " 
             << estDistToSolid << G4endl;
      if( estDistToSolid > 100.0 * motherSolid->GetTolerance() )
      {
        motherSolid->DumpInfo();
        G4Exception("G4NormalNavigation::ComputeStep()",
                    "FarOutsideCurrentVolume", FatalException,
                    "Point is far outside Current Volume !" ); 
      }
      else
        G4Exception("G4NormalNavigation::ComputeStep()", "OutsideCurrentVolume", 
                    JustWarning, "Point is a little outside Current Volume."); 
    }

    // Verification / verbosity
    //
    if ( fVerbose > 1 )
    {
      static G4int precVerf= 20;  // Precision 
      G4int oldprec = G4cout.precision(precVerf);
      G4cout << " G4NormalNavigation::ComputeStep()"
             << " - Information on mother / key daughters ..." << G4endl;
      G4cout << " Type   " << std::setw(12) << "Solid-Name"   << " " 
             << std::setw(3*(6+precVerf))   << " local point" << " "
             << std::setw(4+precVerf)       << "solid-Safety" << " "
             << std::setw(4+precVerf)       << "solid-Step"   << " "
             << std::setw(17)               << "distance Method "
             << std::setw(3*(6+precVerf))   << " local direction" << " "
             << G4endl;
      G4cout << " Mother " << std::setw(12) << motherSolid->GetName() << " "
             << std::setw(4+precVerf)       << localPoint   << " "
             << std::setw(4+precVerf)       << motherSafety << " "
             << G4endl;
      G4cout.precision(oldprec);
    }
}

void
G4NormalNavigation::AlongComputeStepLog(const G4VSolid* sampleSolid,
                                        const G4ThreeVector& samplePoint,
                                        const G4ThreeVector& sampleDirection,
                                        const G4ThreeVector& localDirection,
                                              G4double sampleSafety,
                                              G4double sampleStep,
                                              G4double ourStep)
{
  G4cout << "    Invoked DistanceToIn(p,v) for daughter solid: "
         << sampleSolid->GetName()
         << ". Solid replied: " << sampleStep << G4endl
         << "    For local point p: " << samplePoint << G4endl
         << "    Direction v: " << sampleDirection
         << ", to be considered as 'daughter step'." << G4endl;
  if ( sampleStep<=ourStep )
  {
    // Check to see that the resulting point is indeed in/on volume.
    // This check could eventually be made only for successful candidate.

    if ( sampleStep < kInfinity )
    {
      G4ThreeVector intersectionPoint;
      intersectionPoint= samplePoint + sampleStep * sampleDirection;
      EInside insideIntPt= sampleSolid->Inside(intersectionPoint); 
      G4String solidResponse = "-kInside-";
      if (insideIntPt == kOutside)
        { solidResponse = "-kOutside-"; }
      else if (insideIntPt == kSurface)
        { solidResponse = "-kSurface-"; }
      if( fVerbose == 1 )
      {
        G4cout << "    Invoked Inside() for solid: "
               << sampleSolid->GetName()
               << ". Solid replied: " << solidResponse << G4endl
               << "    For point p: " << intersectionPoint
               << ", considered as 'intersection' point." << G4endl;
      }
      if ( insideIntPt != kSurface )
      {
        G4int oldcoutPrec = G4cout.precision(16); 
        G4cout << "WARNING - G4NormalNavigation::ComputeStep()" << G4endl
               << "          Inaccurate DistanceToIn for solid "
               << sampleSolid->GetName() << G4endl;
        G4cout << "          Solid gave DistanceToIn = " << sampleStep
               << " yet returns " << solidResponse
               << " for this point !" << G4endl; 
        G4cout << "          Point = " << intersectionPoint << G4endl;
        if ( insideIntPt != kInside )
        {
          G4cout << "        DistanceToIn(p) = " 
                 << sampleSolid->DistanceToIn(intersectionPoint)
                 << G4endl;
        }
        if ( insideIntPt != kOutside )
        {
          G4cout << "        DistanceToOut(p) = " 
                 << sampleSolid->DistanceToOut(intersectionPoint)
                 << G4endl;
        }
        G4Exception("G4NormalNavigation::ComputeStep()", 
                    "InaccurateDistanceToIn", JustWarning,
                    "Navigator gets conflicting response from Solid."); 
        G4cout.precision(oldcoutPrec);
      }
    }

    // Verification / verbosity
    //
    if ( fVerbose > 1 )
    {
      static G4int precVerf= 20;  // Precision 
      G4int oldprec = G4cout.precision(precVerf);
      G4cout << " Daught "
             << std::setw(12)         << sampleSolid->GetName() << " "
             << std::setw(4+precVerf) << samplePoint  << " "
             << std::setw(4+precVerf) << sampleSafety << " "
             << std::setw(4+precVerf) << sampleStep   << " "
             << std::setw(16)         << "distanceToIn" << " "
             << std::setw(4+precVerf) << localDirection << " "
             << G4endl;
      G4cout.precision(oldprec);
    }
  }
}

void
G4NormalNavigation::PostComputeStepLog(const G4VSolid* motherSolid,
                                       const G4ThreeVector& localPoint,
                                       const G4ThreeVector& localDirection,
                                             G4double motherStep,
                                             G4double motherSafety)
{
  if ( fCheck )
  {
    if( fVerbose == 1 )
    {
      G4cout << "    Invoked DistanceToOut(p,v,...) for mother solid: "
             << motherSolid->GetName()
             << ". Solid replied: " << motherStep << G4endl
             << "    For local point p: " << localPoint << G4endl
             << "    Direction v: " << localDirection
             << ", to be considered as 'mother step'." << G4endl;
    }
    if( ( motherStep < 0.0 ) || ( motherStep >= kInfinity) )
    {
      G4cerr << "ERROR - G4NormalNavigation::ComputeStep()" << G4endl
             << "        Problem in Navigation"  << G4endl
             << "        Point (local coordinates): "
             << localPoint << G4endl
             << "        Local Direction: " << localDirection << G4endl
             << "        Solid: " << motherSolid->GetName() << G4endl; 
      motherSolid->DumpInfo();
      G4Exception("G4NormalNavigation::ComputeStep()",
                  "PointDistOutInvalid", FatalException,
                  "Current point is outside the current solid !");
    }
  }
  if ( fVerbose > 1 )
  {
    static G4int precVerf= 20;  // Precision 
    G4int oldprec = G4cout.precision(precVerf);
    G4cout << " Mother " << std::setw(12) << motherSolid->GetName() << " "
           << std::setw(4+precVerf)       << localPoint   << " "
           << std::setw(4+precVerf)       << motherSafety << " "
           << std::setw(4+precVerf)       << motherStep   << " "
           << std::setw(16)               << "distanceToOut" << " "
           << std::setw(4+precVerf)       << localDirection << " "
           << G4endl;
    G4cout.precision(oldprec);      
  }
}

void
G4NormalNavigation::ComputeSafetyLog(const G4VSolid* solid,
                                     const G4ThreeVector& point,
                                           G4double safety,
                                           G4bool banner)
{
  G4String volumeType = "daughter";
  if (banner)
  {
    G4cout << "*** G4NormalNavigation::ComputeSafety(): ***" << G4endl;
    volumeType = "mother";
  }
  G4cout << "    Invoked DistanceToOut(p) for "<< volumeType << " solid: "
         << solid->GetName()
         << ". Solid replied: " << safety << G4endl
         << "    For local point p: " << point
         << ", to be considered as '"<< volumeType << " safety'." << G4endl;
}
