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
// class G4NavigationLogger Implementation
//
// Author: G.Cosmo, 2010
// --------------------------------------------------------------------

#include <iomanip>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4NavigationLogger.hh"
#include "G4GeometryTolerance.hh"

using CLHEP::millimeter;

G4NavigationLogger::G4NavigationLogger(const G4String& id)
  : fId(id)
{
}

G4NavigationLogger::~G4NavigationLogger()
{
}

// ********************************************************************
// PreComputeStepLog
// ********************************************************************
//
void
G4NavigationLogger::PreComputeStepLog(const G4VPhysicalVolume* motherPhysical,
                                            G4double motherSafety,
                                      const G4ThreeVector& localPoint) const
{
  G4VSolid* motherSolid = motherPhysical->GetLogicalVolume()->GetSolid();
  G4String fType = fId + "::ComputeStep()";
  
  if ( fVerbose == 1 || fVerbose > 4 )       
  {
    G4cout << "*************** " << fType << " *****************" << G4endl
           << " VolType "
           << std::setw(15) << "Safety/mm" << " "
           << std::setw(15) << "Distance/mm" << " "
           << std::setw(52) << "Position (local coordinates)"
           << " - Solid" << G4endl;
    G4cout << "  Mother "
           << std::setw(15) << motherSafety / millimeter << " " 
           << std::setw(15) << "N/C"        << " " << localPoint << " - "
           << motherSolid->GetEntityType() << ": " << motherSolid->GetName()
           << G4endl;
  }
  if ( motherSafety < 0.0 )
  {
    std::ostringstream message;
    message << "Negative Safety In Voxel Navigation !" << G4endl
            << "        Current solid " << motherSolid->GetName()
            << " gave negative safety: " << motherSafety / millimeter << G4endl
            << "        for the current (local) point " << localPoint;
    message << " Solid info: " << *motherSolid << G4endl;      
    G4Exception(fType, "GeomNav0003", FatalException, message);
  }
  if( motherSolid->Inside(localPoint)==kOutside )
  {
    std::ostringstream message;
    message << "Point is outside Current Volume - " << G4endl
            << "          Point " << localPoint / millimeter
            << " is outside current volume '" << motherPhysical->GetName()
            << "'" << G4endl;
    G4double estDistToSolid= motherSolid->DistanceToIn(localPoint); 
    message << "          Estimated isotropic distance to solid (distToIn)= " 
            << estDistToSolid << G4endl;
    if( estDistToSolid > 100.0 * motherSolid->GetTolerance() )
    {
      message << " Solid info: " << *motherSolid << G4endl;
      G4Exception(fType, "GeomNav0003", JustWarning, message,
                  "Point is far outside Current Volume !" ); 
    }
    else
    {
      G4Exception(fType, "GeomNav1001", JustWarning, message,
                  "Point is a little outside Current Volume.");
    }
  }

  // Verification / verbosity
  //
  if ( fVerbose > 1 )
  {
    static const G4int precVerf = 16;  // Precision 
    G4long oldprec = G4cout.precision(precVerf);
    G4cout << " - Information on mother / key daughters ..." << G4endl;
    G4cout << "  Type   " << std::setw(12) << "Solid-Name"   << " " 
           << std::setw(3*(6+precVerf))    << " local point" << " "
           << std::setw(4+precVerf)        << "solid-Safety" << " "
           << std::setw(4+precVerf)        << "solid-Step"   << " "
           << std::setw(17)                << "distance Method "
           << std::setw(3*(6+precVerf))    << " local direction" << " "
           << G4endl;
    G4cout << "  Mother " << std::setw(12) << motherSolid->GetName() << " "
           << std::setw(4+precVerf)        << localPoint   << " "
           << std::setw(4+precVerf)        << motherSafety << " "
           << G4endl;
    G4cout.precision(oldprec);
  }
}

// ********************************************************************
// AlongComputeStepLog
// ********************************************************************
//
void
G4NavigationLogger::AlongComputeStepLog(const G4VSolid* sampleSolid,
                                        const G4ThreeVector& samplePoint,
                                        const G4ThreeVector& sampleDirection,
                                        const G4ThreeVector& localDirection,
                                              G4double sampleSafety,
                                              G4double sampleStep) const
{
  // Check to see that the resulting point is indeed in/on volume.
  // This check could eventually be made only for successful candidate.

  if ( sampleStep < kInfinity )
  {
    G4ThreeVector intersectionPoint;
    intersectionPoint = samplePoint + sampleStep * sampleDirection;
    EInside insideIntPt = sampleSolid->Inside(intersectionPoint); 
    G4String fType = fId + "::ComputeStep()";

    G4String solidResponse = "-kInside-";
    if (insideIntPt == kOutside)
      { solidResponse = "-kOutside-"; }
    else if (insideIntPt == kSurface)
      { solidResponse = "-kSurface-"; }

    if ( fVerbose == 1 || fVerbose > 4 )
    {
      G4cout << "    Invoked Inside() for solid: "
             << sampleSolid->GetName()
             << ". Solid replied: " << solidResponse << G4endl
             << "    For point p: " << intersectionPoint
             << ", considered as 'intersection' point." << G4endl;
    }

    G4double safetyIn = -1, safetyOut = -1;  //  Set to invalid values
    G4double newDistIn = -1,  newDistOut = -1;
    if( insideIntPt != kInside )
    {
      safetyIn = sampleSolid->DistanceToIn(intersectionPoint);
      newDistIn = sampleSolid->DistanceToIn(intersectionPoint,
                                            sampleDirection);
    }
    if( insideIntPt != kOutside )
    {
      safetyOut = sampleSolid->DistanceToOut(intersectionPoint);
      newDistOut = sampleSolid->DistanceToOut(intersectionPoint,
                                              sampleDirection);
    }
    if( insideIntPt != kSurface )
    {
      std::ostringstream message;
      message.precision(16); 
      message << "Conflicting response from Solid." << G4endl
              << "          Inaccurate solid DistanceToIn"
              << " for solid " << sampleSolid->GetName() << G4endl
              << "          Solid gave DistanceToIn = "
              << sampleStep << " yet returns " << solidResponse
              << " for this point !" << G4endl
              << "          Original Point     = " << samplePoint << G4endl
              << "          Original Direction = " << sampleDirection << G4endl              
              << "          Intersection Point = " << intersectionPoint << G4endl
              << "            Safety values: " << G4endl;
      if ( insideIntPt != kInside )
      {
        message << "          DistanceToIn(p)  = " << safetyIn;
      }
      if ( insideIntPt != kOutside )
      {
        message << "          DistanceToOut(p) = " << safetyOut;
      }
      message << G4endl;
      message << " Solid Parameters: " << *sampleSolid;
      G4Exception(fType, "GeomNav1001", JustWarning, message);
    }
    else
    {  
      // If it is on the surface, *ensure* that either DistanceToIn
      // or DistanceToOut returns a finite value ( >= Tolerance).
      //
      if( std::max( newDistIn, newDistOut ) <=
          G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() )
      { 
        std::ostringstream message;
        message << "Zero from both Solid DistanceIn and Out(p,v)." << G4endl
                << "  Identified point for which the solid " 
                << sampleSolid->GetName() << G4endl
                << "  has MAJOR problem:  " << G4endl
                << "  --> Both DistanceToIn(p,v) and DistanceToOut(p,v) "
                << "return Zero, an equivalent value or negative value."
                << G4endl 
                << "    Solid: " << sampleSolid << G4endl
                << "    Point p= " << intersectionPoint << G4endl
                << "    Direction v= " << sampleDirection << G4endl
                << "    DistanceToIn(p,v)     = " << newDistIn << G4endl
                << "    DistanceToOut(p,v,..) = " << newDistOut << G4endl
                << "    Safety values: " << G4endl
                << "      DistanceToIn(p)  = " << safetyIn << G4endl
                << "      DistanceToOut(p) = " << safetyOut;
        G4Exception(fType, "GeomNav0003", FatalException, message);
      }
    }

    // Verification / verbosity
    //
    if ( fVerbose > 1 )
    {
      static const G4int precVerf= 20;  // Precision 
      G4long oldprec = G4cout.precision(precVerf);
      G4cout << "Daughter "
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

// ********************************************************************
// CheckDaughterEntryPoint
// ********************************************************************
//
void
G4NavigationLogger::CheckDaughterEntryPoint(const G4VSolid* sampleSolid,
                                            const G4ThreeVector& samplePoint,
                                            const G4ThreeVector& sampleDirection,
                                            const G4VSolid* motherSolid,
                                            const G4ThreeVector& localPoint,
                                            const G4ThreeVector& localDirection,
                                                  G4double motherStep,
                                                  G4double sampleStep) const
{
  const G4double kCarTolerance = motherSolid->GetTolerance();
   
  // Double check the expected condition of being called
  //
  G4bool SuspiciousDaughterDist = ( sampleStep >= motherStep )
                               && ( sampleStep < kInfinity );

  if( sampleStep >= kInfinity )
  {
    G4ExceptionDescription msg;
    msg.precision(12);
    msg << " WARNING - Called with 'infinite' step. " << G4endl;
    msg << "    Checks have no meaning if daughter step is infinite." << G4endl;
    msg << "    kInfinity  = " << kInfinity  / millimeter << G4endl;
    msg << "    sampleStep = " << sampleStep / millimeter << G4endl;
    msg << "    sampleStep < kInfinity " << (sampleStep<kInfinity) << G4endl;
    msg << "    kInfinity - sampleStep " << (kInfinity-sampleStep) / millimeter << G4endl;
    msg << " Returning immediately.";
    G4Exception("G4NavigationLogger::CheckDaughterEntryPoint()",
                "GeomNav0003", JustWarning, msg);
    return; 
  }

  // The intersection point with the daughter is after the exit point
  // from the mother volume !!
  // This is legal / allowed to occur only if the mother is concave
  // ****************************************************************
  // If mother is convex the daughter volume must be encountered
  // before the exit from the current volume!
   
  // Check #1) whether the track will re-enter the current mother 
  //           in the extension past its current exit point 
  G4ThreeVector localExitMotherPos = localPoint+motherStep*localDirection;
  G4double distExitToReEntry = motherSolid->DistanceToIn(localExitMotherPos,
                                                         localDirection); 
   
  // Check #2) whether the 'entry' point in the daughter is inside the mother
  //
  G4ThreeVector localEntryInDaughter = localPoint+sampleStep*localDirection; 
  EInside insideMother = motherSolid->Inside( localEntryInDaughter ); 
   
  G4String solidResponse = "-kInside-";
  if (insideMother == kOutside)       { solidResponse = "-kOutside-"; }
  else if (insideMother == kSurface)  { solidResponse = "-kSurface-"; }

  G4double       distToReEntry = distExitToReEntry + motherStep;
  G4ThreeVector  localReEntryPoint = localPoint+distToReEntry*localDirection;

  // Clear error  -- Daughter entry point is bad
  constexpr G4double eps= 1.0e-10;
  G4bool DaughterEntryIsOutside = SuspiciousDaughterDist
     && ( (sampleStep * (1.0+eps) < distToReEntry) || (insideMother == kOutside ) );
  G4bool EntryIsMotherExit = std::fabs(sampleStep-motherStep) < kCarTolerance;

  // Check for more subtle error - is exit point of daughter correct ?
  G4ThreeVector sampleEntryPoint = samplePoint+sampleStep*sampleDirection;
  G4double sampleCrossingDist = sampleSolid->DistanceToOut( sampleEntryPoint,
                                                            sampleDirection );
  G4double      sampleExitDist = sampleStep+sampleCrossingDist;
  G4ThreeVector sampleExitPoint = samplePoint+sampleExitDist*sampleDirection;

  G4bool TransitProblem = ( (sampleStep < motherStep)
                         && (sampleExitDist > motherStep + kCarTolerance) )
           || ( EntryIsMotherExit && (sampleCrossingDist > kCarTolerance) );

  if( DaughterEntryIsOutside
       || TransitProblem
       || (SuspiciousDaughterDist && (fVerbose > 3) ) )
  {
    G4ExceptionDescription msg;
    msg.precision(16);

    if( DaughterEntryIsOutside )
    {   
      msg << "WARNING> Intersection distance to Daughter volume is further"
          <<    " than the distance to boundary." << G4endl
          << "  It appears that part of the daughter volume is *outside*"
          <<    " this mother. " << G4endl;
      msg << "  One of the following checks signaled a problem:" << G4endl
          << "  -sampleStep (dist to daugh) <  mother-exit dist + distance "
          <<      "to ReEntry point for mother " << G4endl
          << "  -position of daughter intersection is outside mother volume."
          << G4endl;
    }
    else if( TransitProblem )
    {
      G4double protrusion = sampleExitDist - motherStep;

      msg << "WARNING>  Daughter volume extends beyond mother boundary. "
          << G4endl;
      if ( ( sampleStep < motherStep )
        && (sampleExitDist > motherStep + kCarTolerance ) )
      {
        // 1st Issue with Daughter
        msg << "        Crossing distance in the daughter causes is to extend"
            << " beyond the mother exit. " << G4endl;
        msg << "        Length protruding = " << protrusion << G4endl;
      }
      if( EntryIsMotherExit )
      {
        // 1st Issue with Daughter
        msg << "        Intersection distance to Daughter is within "
            << " tolerance of the distance" << G4endl;
        msg << "        to the mother boundary * and * " << G4endl;
        msg << "        the crossing distance in the daughter is > tolerance."
            << G4endl;           
      }         
    }      
    else
    {
      msg << "NearMiss> Intersection to Daughter volume is in extension past the"
          <<   " current exit point of the mother volume." << G4endl;
      msg << "          This is not an error - just an unusual occurrence,"
          <<   " possible in the case of concave volume. " << G4endl;
    }
    msg << "---- Information about intersection with daughter, mother: "
        << G4endl;
    msg << "    sampleStep (daughter) = " << sampleStep << G4endl
        << "    motherStep            = " << motherStep << G4endl
        << "    distToRentry(mother)  = " << distToReEntry << G4endl
        << "    Inside(entry pnt daug): " << solidResponse << G4endl
        << "    dist across daughter  = " << sampleCrossingDist << G4endl;
    msg << " Mother Name (Solid) : " << motherSolid->GetName() << G4endl
        << " In local (mother) coordinates: " << G4endl
        << "    Starting     Point    = " << localPoint << G4endl
        << "    Direction             = " << localDirection << G4endl
        << "    Exit Point    (mother)= " << localExitMotherPos << G4endl
        << "    Entry Point (daughter)= " << localPoint+sampleStep*localDirection
        << G4endl;
    if( distToReEntry < kInfinity )
    {
      msg << "    ReEntry Point (mother)= " << localReEntryPoint << G4endl;
    }
    else
    {
      msg << "    No ReEntry - track does not encounter mother volume again! "
          << G4endl;
    }
    msg << " Daughter Name (Solid): " << sampleSolid->GetName() << G4endl
        << " In daughter coordinates: " << G4endl
        << "    Starting     Point    = " << samplePoint << G4endl
        << "    Direction             = " << sampleDirection << G4endl
        << "    Entry Point (daughter)= " << sampleEntryPoint
        << G4endl;
    msg << "  Description of mother solid: " << G4endl
        << *motherSolid << G4endl
        << "  Description of daughter solid: " << G4endl
        << *sampleSolid << G4endl;
    G4String fType = fId + "::ComputeStep()";

    if( DaughterEntryIsOutside || TransitProblem )
    {
      G4Exception(fType, "GeomNav0003", JustWarning, msg);
    }
    else
    {
      G4cout << fType
             << " -- Checked distance of Entry to daughter vs exit of mother"
             << G4endl;
      G4cout << msg.str();
      G4cout << G4endl;
    }
  }
}

// ********************************************************************
// PostComputeStepLog
// ********************************************************************
//
void
G4NavigationLogger::PostComputeStepLog(const G4VSolid* motherSolid,
                                       const G4ThreeVector& localPoint,
                                       const G4ThreeVector& localDirection,
                                             G4double motherStep,
                                             G4double motherSafety) const
{
  if ( fVerbose == 1 || fVerbose > 4 )     
  {
    G4cout << "  Mother "
           << std::setw(15) << motherSafety << " " 
           << std::setw(15) << motherStep   << " " << localPoint   << " - "
           << motherSolid->GetEntityType() << ": " << motherSolid->GetName()
           << G4endl;
  }
  if( ( motherStep < 0.0 ) || ( motherStep >= kInfinity) )
  {
    G4String fType = fId + "::ComputeStep()";
    G4long oldPrOut = G4cout.precision(16); 
    G4long oldPrErr = G4cerr.precision(16);
    std::ostringstream message;
    message << "Current point is outside the current solid !" << G4endl
            << "        Problem in Navigation"  << G4endl
            << "        Point (local coordinates): "
            << localPoint << G4endl
            << "        Local Direction: " << localDirection << G4endl
            << "        Solid: " << motherSolid->GetName(); 
    motherSolid->DumpInfo();
    G4Exception(fType, "GeomNav0003", FatalException, message);
    G4cout.precision(oldPrOut);
    G4cerr.precision(oldPrErr);
  }
  if ( fVerbose > 1 )
  {
    static const G4int precVerf = 20;  // Precision 
    G4long oldprec = G4cout.precision(precVerf);
    G4cout << "  Mother " << std::setw(12) << motherSolid->GetName() << " "
           << std::setw(4+precVerf)       << localPoint   << " "
           << std::setw(4+precVerf)       << motherSafety << " "
           << std::setw(4+precVerf)       << motherStep   << " "
           << std::setw(16)               << "distanceToOut" << " "
           << std::setw(4+precVerf)       << localDirection << " "
           << G4endl;
    G4cout.precision(oldprec);      
  }
}

// ********************************************************************
// ComputeSafetyLog
// ********************************************************************
//
void
G4NavigationLogger::ComputeSafetyLog(const G4VSolid* solid,
                                     const G4ThreeVector& point,
                                           G4double safety,
                                           G4bool isMotherVolume, 
                                           G4int banner) const
{
  if( banner < 0 )
  {
    banner = isMotherVolume;
  }
  if( fVerbose >= 1 )
  {
    G4String volumeType = isMotherVolume ? " Mother " : "Daughter";
    if (banner)
    {
      G4cout << "************** " << fId << "::ComputeSafety() ****************"
             << G4endl;
      G4cout << " VolType "
             << std::setw(15) << "Safety/mm" << " "
             << std::setw(52) << "Position (local coordinates)"
             << " - Solid" << G4endl;
    }
    G4cout << volumeType
           << std::setw(15) << safety << " " << point  << " - "
           << solid->GetEntityType() << ": " << solid->GetName() << G4endl;
  }
}

// ********************************************************************
// PrintDaughterLog
// ********************************************************************
//
void
G4NavigationLogger::PrintDaughterLog (const G4VSolid* sampleSolid,
                                      const G4ThreeVector& samplePoint,
                                            G4double sampleSafety,
                                            G4bool   withStep,
                                      const G4ThreeVector& sampleDirection,                                                                               G4double sampleStep ) const
{
  if ( fVerbose >= 1 )
  {
    G4long oldPrec = G4cout.precision(8);
    G4cout << "Daughter "
           << std::setw(15) << sampleSafety << " ";
    if (withStep)  // (sampleStep != -1.0 )
    {
      G4cout << std::setw(15) << sampleStep << " ";
    }
    else
    {
      G4cout << std::setw(15) << "Not-Available" << " ";
    }
    G4cout << samplePoint   << " - "
           << sampleSolid->GetEntityType() << ": " << sampleSolid->GetName();
    if( withStep )
    {
      G4cout << " dir= " << sampleDirection;
    }
    G4cout << G4endl;
    G4cout.precision(oldPrec);
  }
}

// ********************************************************************
// CheckAndReportBadNormal
// ********************************************************************
//
G4bool
G4NavigationLogger::
CheckAndReportBadNormal(const G4ThreeVector& unitNormal,
                        const G4ThreeVector& localPoint,
                        const G4ThreeVector& localDirection,
                              G4double         step,
                        const G4VSolid*      solid,
                        const char*          msg ) const
{
  G4double normMag2 = unitNormal.mag2();
  G4bool badLength = ( std::fabs ( normMag2 - 1.0 ) > CLHEP::perMillion );

  if( badLength )
  {
    G4double  normMag = std::sqrt(normMag2); 
    G4ExceptionDescription message;
    message.precision(10);
    message << "============================================================"
            << G4endl;
    message << " WARNING>  Normal is not a unit vector. "
            << "  - but |normal|   = "  << normMag
            << "  - and |normal|^2     = "  << normMag2 << G4endl
            << "    which differ from 1.0 by: " <<  G4endl 
            << "        |normal|-1 = " << normMag - 1.0
            << "    and |normal|^2 - 1 = " << normMag2 - 1.0 << G4endl
            << "   n = " << unitNormal << G4endl;
    message << " Info string: " << msg << G4endl;
    message << "============================================================"
            << G4endl;

    message.precision(16);

    message << " Information on call to DistanceToOut: " << G4endl;
    message << "   Position  = " << localPoint << G4endl
            << "   Direction = " << localDirection << G4endl;
    message << "   Obtained> distance      = " << step << G4endl;
    message << "           > Exit position = " << localPoint+step*localDirection
            << G4endl;
    message << " Parameters of solid:     " << G4endl;
    message << *solid; 
    message << "============================================================";
      
    G4String fMethod = fId + "::ComputeStep()";
    G4Exception( fMethod, "GeomNav0003", JustWarning, message); 
  }
  return badLength;
}

// ********************************************************************
// CheckAndReportBadNormal - due to Rotation Matrix
// ********************************************************************
//
G4bool
G4NavigationLogger::
CheckAndReportBadNormal(const G4ThreeVector& rotatedNormal,
                        const G4ThreeVector& originalNormal,
                        const G4RotationMatrix& rotationM,
                        const char*          msg ) const
{
  G4double  normMag2 = rotatedNormal.mag2();
  G4bool badLength = ( std::fabs ( normMag2 - 1.0 ) > CLHEP::perMillion );

  if( badLength )
  {
    G4double  normMag = std::sqrt(normMag2); 
    G4ExceptionDescription message;
    message.precision(10);
    message << "============================================================"
            << G4endl;
    message << " WARNING>  Rotated n(ormal) is not a unit vector. " << G4endl
            << "     |normal|   = "  << normMag
            << "   and |normal|^2     = "  << normMag2 << G4endl
            << "   Diff from 1.0: " <<  G4endl 
            << "     |normal|-1 = " << normMag - 1.0
            << "   and |normal|^2 - 1 = " << normMag2 - 1.0 << G4endl;
    message << "   Rotated  n = (" << rotatedNormal.x() << "," << rotatedNormal.y() << "," 
            << rotatedNormal.z() << ")" << G4endl;
    message << "   Original n = (" << originalNormal.x() << "," << originalNormal.y() << "," 
            << originalNormal.z() << ")" << G4endl;
    message << " Info string: " << msg << G4endl;
    message << "============================================================"
            << G4endl;

    message.precision(16);

    message << " Information on RotationMatrix : " << G4endl;
    message << " Original: " << G4endl;
    message << rotationM << G4endl;
    message << " Inverse (used in transformation): " << G4endl;
    message << rotationM.inverse() << G4endl;    
    message << "============================================================";
      
    G4String fMethod = fId + "::ComputeStep()";
    G4Exception( fMethod, "GeomNav0003", JustWarning, message); 
  }
  return badLength;
}

// ********************************************************************
// ReportOutsideMother
// ********************************************************************
//
// Report Exception if point is outside mother.
// Fatal exception will be used if either 'checkMode error is > triggerDist
//
void
G4NavigationLogger::ReportOutsideMother(const G4ThreeVector& localPoint,
                                        const G4ThreeVector& localDirection,
                                        const G4VPhysicalVolume* physical,
                                              G4double triggerDist) const                                   
{
  const G4LogicalVolume* logicalVol = physical
                                    ? physical->GetLogicalVolume() : nullptr;
  const G4VSolid* solid = logicalVol
                        ? logicalVol->GetSolid() : nullptr;

  G4String fMethod = fId + "::ComputeStep()";

  if( solid == nullptr )
  {
    G4Exception(fMethod, "GeomNav0003", FatalException,
                "Erroneous call to ReportOutsideMother: no Solid is available");
    return;
  }
  const G4double kCarTolerance = solid->GetTolerance();

  // Double check reply - it should be kInfinity
  const G4double distToOut = solid->DistanceToOut(localPoint, localDirection);  
  const EInside  inSolid   = solid->Inside(localPoint);
  const G4double safetyToIn  = solid->DistanceToIn(localPoint);
  const G4double safetyToOut = solid->DistanceToOut(localPoint);
  // const G4double distToInPos =
  //                solid->DistanceToIn(localPoint, localDirection);

  // 1. Check consistency between Safety obtained and report from distance
  //     We must ensure that (mother)Safety <= 0.0
  //       in the case that the point is outside!
  //    [ If this is not the case, this problem can easily go undetected,
  //       except in Check mode ! ] 
  if( safetyToOut > kCarTolerance
      && ( distToOut < 0.0 || distToOut >= kInfinity ) )
  {
     G4ExceptionDescription msg1;     
     // fNavClerk->ReportBadSafety(localPoint, localDirection,
     //                     motherPhysical, motherSafety, motherStep);
     msg1 << " Dangerous inconsistency in response of solid." << G4endl
          << "    Solid type: " << solid->GetEntityType()
          << "    Name= " << solid->GetName()           << G4endl;
     msg1 << " Mother volume gives safety > 0 despite being called for *Outside* point "
          << G4endl
          << "   Location = " << localPoint     << G4endl
          << "   Direction= " << localDirection << G4endl
          << "   - Safety (Isotropic d) = " << safetyToOut   << G4endl
          << "   - Intersection Distance= " << distToOut << G4endl
          << G4endl;
     G4Exception( fMethod, "GeomNav0123", JustWarning, msg1);
  }

  // 2. Inconsistency - Too many distances are zero (or will be rounded to zero)

//  if( std::fabs(distToOut) < kCarTolerance
//   && std::fabs(distToInPos) < kCarTolerance ) 
//  {
     // If both distanceToIn and distanceToOut (p,v) are zero for
     // one direction, the particle could get stuck!
//  }

  G4ExceptionDescription msg;
  msg.precision(10);
  
  if( std::fabs(distToOut) < kCarTolerance )
  {
    // 3. Soft error - safety is not rounded to zero
    //    Report nothing - except in 'loud' mode
    if( fReportSoftWarnings )
    {
      msg << " Warning>  DistanceToOut(p,v): "
          << "Distance from surface is not rounded to zero" << G4endl;
    }
    else
    { 
      return;
    }
  }
  else
  {
    // 4. General message - complain that the point is outside!
    //     and provide all information about the current location,
    //     direction and the answers of the solid
    msg << "============================================================"
        << G4endl;
    msg << " WARNING>  Current Point appears to be Outside mother volume !! "
        << G4endl;
    msg << "   Response of DistanceToOut was negative or kInfinity"
        << " when called in " << fMethod << G4endl;
  }

  // Generate and 'print'/stream all the information needed
  this->ReportVolumeAndIntersection(msg, localPoint, localDirection, physical); 
  
  // Default for distance of 'major' error
  if( triggerDist <= 0.0 )
  {
    triggerDist = std::max ( 1.0e+6 * kCarTolerance,  // Well beyond tolerance
                             fMinTriggerDistance );
  }

  G4bool majorError = inSolid == kOutside
                    ? ( safetyToIn > triggerDist )
                    : ( safetyToOut > triggerDist );
  
  G4ExceptionSeverity exceptionType = JustWarning;
  if ( majorError )
  {
    exceptionType = FatalException;
  }
  
  G4Exception( fMethod, "GeomNav0003", exceptionType, msg);
}

namespace  G4NavigationLogger_Namespace
{
  const G4String EInsideNames[3] = { "kOutside", "kSurface", "kInside" }; 
}

void G4NavigationLogger::
ReportVolumeAndIntersection( std::ostream& os,
                             const G4ThreeVector&     localPoint,
                             const G4ThreeVector&     localDirection,
                             const G4VPhysicalVolume* physical ) const 
{
  G4String fMethod = fId + "::ComputeStep()";   
  const G4LogicalVolume* logicalVol = physical
                                    ? physical->GetLogicalVolume() : nullptr;
  const G4VSolid* solid = logicalVol
                        ? logicalVol->GetSolid() : nullptr;
  if( solid == nullptr )
  {
     os << " ERROR> Solid is not available. Logical Volume = "
        << logicalVol << std::endl;
     return;
  }
  const G4double kCarTolerance = solid->GetTolerance();

  // Double check reply - it should be kInfinity
  const G4double distToOut = solid->DistanceToOut(localPoint, localDirection);
  const G4double distToOutNeg = solid->DistanceToOut(localPoint,
                                                    -localDirection);
  const EInside  inSolid   = solid->Inside(localPoint);
  const G4double safetyToIn  = solid->DistanceToIn(localPoint);
  const G4double safetyToOut = solid->DistanceToOut(localPoint);

  const G4double distToInPos = solid->DistanceToIn(localPoint, localDirection);
  const G4double distToInNeg = solid->DistanceToIn(localPoint, -localDirection);  

  const G4ThreeVector exitNormal = solid->SurfaceNormal(localPoint); 

  // Double check whether points nearby are in/surface/out
  const G4double epsilonDist = 1000.0 * kCarTolerance;
  const G4ThreeVector PointPlusDir = localPoint + epsilonDist * localDirection;
  const G4ThreeVector PointMinusDir = localPoint - epsilonDist * localDirection;
  const G4ThreeVector PointPlusNorm = localPoint + epsilonDist * exitNormal;  
  const G4ThreeVector PointMinusNorm = localPoint - epsilonDist * exitNormal;

  const EInside inPlusDir = solid->Inside(PointPlusDir);
  const EInside inMinusDir = solid->Inside(PointMinusDir);  
  const EInside inPlusNorm = solid->Inside(PointPlusNorm);
  const EInside inMinusNorm = solid->Inside(PointMinusNorm);  

  // Basic information 
  os << "   Current physical volume = " << physical->GetName() << G4endl;
  os << "   Position (loc)  = " << localPoint << G4endl
     << "   Direction (dir) = " << localDirection << G4endl;
  os << " For confirmation:" << G4endl;
  os << "   Response of DistanceToOut (loc, +dir)= " << distToOut << G4endl;
  os << "   Response of DistanceToOut (loc, -dir)= " << distToOutNeg << G4endl;
  
  os << "   Inside responds = " << inSolid << " , ie: ";
  if( inSolid == kOutside )
  {
    os << " Outside -- a problem, as observed in " << fMethod << G4endl;
  }
  else if( inSolid == kSurface )
  {
    os << " Surface -- unexpected / inconsistent response ! " << G4endl;
  }
  else
  {
    os << " Inside  -- unexpected / inconsistent response ! " << G4endl;
  }
  os << "   Obtain safety(ToIn) = " << safetyToIn << G4endl;
  os << "   Obtain safety(ToOut) = " << safetyToOut << G4endl;
  os << " Response of DistanceToIn (loc, +dir)= "  << distToInPos << G4endl;
  os << " Response of DistanceToIn (loc, -dir)= "  << distToInNeg << G4endl;

  os << " Exit Normal at loc = " << exitNormal << G4endl;
  os << "     Dir . Normal   = " << exitNormal.dot( localDirection );
  os << G4endl;

  os << " Checking points moved from position by distance/direction." << G4endl
     << " Solid responses: " << G4endl
     << "  +eps in direction :    "
     << G4NavigationLogger_Namespace::EInsideNames[inPlusDir] 
     << "  +eps in Normal  :    "
     << G4NavigationLogger_Namespace::EInsideNames[inPlusNorm]  << G4endl
     << "  -eps in direction :    "
     << G4NavigationLogger_Namespace::EInsideNames[inMinusDir]
     << "  -eps in Normal  :    "
     << G4NavigationLogger_Namespace::EInsideNames[inMinusNorm]  << G4endl;
     
  os << " Parameters of solid:     " << G4endl;
  os << *solid;
  os << "============================================================";
}
