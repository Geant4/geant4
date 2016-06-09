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
// $Id$
//
//
// class G4NavigationLogger Implementation
//
// Author: G.Cosmo, 2010
//
// --------------------------------------------------------------------

#include <iomanip>

#include "G4NavigationLogger.hh"
#include "G4GeometryTolerance.hh"

G4NavigationLogger::G4NavigationLogger(const G4String& id)
  : fId(id), fVerbose(0)
{
}

G4NavigationLogger::~G4NavigationLogger()
{
}

void
G4NavigationLogger::PreComputeStepLog(const G4VPhysicalVolume* motherPhysical,
                                            G4double motherSafety,
                                      const G4ThreeVector& localPoint) const
{
    G4VSolid* motherSolid = motherPhysical->GetLogicalVolume()->GetSolid();
    G4String fType = fId + "::ComputeStep()";
    if( fVerbose == 1 )
    {
      G4cout << "*************** " << fType << " *****************" << G4endl
             << " VolType "
             << std::setw(15) << "Safety/mm" << " "
             << std::setw(15) << "Distance/mm" << " "
             << std::setw(52) << "Position (local coordinates)"
             << " - Solid" << G4endl;
      G4cout << "  Mother "
             << std::setw(15) << motherSafety << " " 
             << std::setw(15) << "N/C"        << " " << localPoint << " - "
             << motherSolid->GetEntityType() << ": " << motherSolid->GetName()
             << G4endl;
    }
    if ( motherSafety < 0.0 )
    {
      std::ostringstream message;
      message << "Negative Safety In Voxel Navigation !" << G4endl
             << "        Current solid " << motherSolid->GetName()
             << " gave negative safety: " << motherSafety << G4endl
             << "        for the current (local) point " << localPoint;
      motherSolid->DumpInfo();
      G4Exception(fType, "GeomNav0003", FatalException, message);
    }
    if( motherSolid->Inside(localPoint)==kOutside )
    {
      std::ostringstream message;
      message << "Point is outside Current Volume - " << G4endl
              << "          Point " << localPoint
              << " is outside current volume " << motherPhysical->GetName()
              << G4endl;
      G4double estDistToSolid= motherSolid->DistanceToIn(localPoint); 
      message << "          Estimated isotropic distance to solid (distToIn)= " 
              << estDistToSolid;
      if( estDistToSolid > 100.0 * motherSolid->GetTolerance() )
      {
        motherSolid->DumpInfo();
        G4Exception(fType, "GeomNav0003", FatalException, message,
                    "Point is far outside Current Volume !" ); 
      }
      else
        G4Exception(fType, "GeomNav1001", JustWarning, message,
                    "Point is a little outside Current Volume."); 
    }

    // Verification / verbosity
    //
    if ( fVerbose > 1 )
    {
      static G4int precVerf= 20;  // Precision 
      G4int oldprec = G4cout.precision(precVerf);
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
    intersectionPoint= samplePoint + sampleStep * sampleDirection;
    EInside insideIntPt= sampleSolid->Inside(intersectionPoint); 
    G4String fType = fId + "::ComputeStep()";

    G4String solidResponse = "-kInside-";
    if (insideIntPt == kOutside)
      { solidResponse = "-kOutside-"; }
    else if (insideIntPt == kSurface)
      { solidResponse = "-kSurface-"; }

    if ( fVerbose == 1 )
    {
      G4cout << "    Invoked Inside() for solid: "
             << sampleSolid->GetName()
             << ". Solid replied: " << solidResponse << G4endl
             << "    For point p: " << intersectionPoint
             << ", considered as 'intersection' point." << G4endl;
    }

    G4double safetyIn= -1, safetyOut= -1;  //  Set to invalid values
    G4double newDistIn= -1,  newDistOut= -1;
    if( insideIntPt != kInside )
    {
      safetyIn= sampleSolid->DistanceToIn(intersectionPoint);
      newDistIn= sampleSolid->DistanceToIn(intersectionPoint,
                                           sampleDirection);
    }
    if( insideIntPt != kOutside )
    {
      safetyOut= sampleSolid->DistanceToOut(intersectionPoint);
      newDistOut= sampleSolid->DistanceToOut(intersectionPoint,
                                             sampleDirection);
    }
    if( insideIntPt != kSurface )
    {
      G4int oldcoutPrec = G4cout.precision(16); 
      std::ostringstream message;
      message << "Conflicting response from Solid." << G4endl
              << "          Inaccurate solid DistanceToIn"
              << " for solid " << sampleSolid->GetName() << G4endl
              << "          Solid gave DistanceToIn = "
              << sampleStep << " yet returns " << solidResponse
              << " for this point !" << G4endl
              << "          Point = " << intersectionPoint << G4endl
              << "          Safety values: " << G4endl;
      if ( insideIntPt != kInside )
      {
        message << "          DistanceToIn(p)  = " << safetyIn;
      }
      if ( insideIntPt != kOutside )
      {
        message << "          DistanceToOut(p) = " << safetyOut;
      }
      G4Exception(fType, "GeomNav1001", JustWarning, message);
      G4cout.precision(oldcoutPrec);
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
      static G4int precVerf= 20;  // Precision 
      G4int oldprec = G4cout.precision(precVerf);
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

void
G4NavigationLogger::PostComputeStepLog(const G4VSolid* motherSolid,
                                       const G4ThreeVector& localPoint,
                                       const G4ThreeVector& localDirection,
                                             G4double motherStep,
                                             G4double motherSafety) const
{
  if( fVerbose == 1 )
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
    G4int oldPrOut= G4cout.precision(16); 
    G4int oldPrErr= G4cerr.precision(16);
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
    static G4int precVerf= 20;  // Precision 
    G4int oldprec = G4cout.precision(precVerf);
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

void
G4NavigationLogger::ComputeSafetyLog(const G4VSolid* solid,
                                     const G4ThreeVector& point,
                                           G4double safety,
                                           G4bool banner) const
{
  G4String volumeType = "Daughter ";
  if (banner)
  {
    G4cout << "************** " << fId << "::ComputeSafety() ****************" << G4endl;
    G4cout << " VolType "
           << std::setw(15) << "Safety/mm" << " "
           << std::setw(52) << "Position (local coordinates)"
           << " - Solid" << G4endl;
    volumeType = "  Mother ";
  }
  G4cout << volumeType
         << std::setw(15) << safety << " " << point  << " - "
         << solid->GetEntityType() << ": " << solid->GetName() << G4endl;
}

void
G4NavigationLogger::PrintDaughterLog (const G4VSolid* sampleSolid,
                                      const G4ThreeVector& samplePoint,
                                            G4double sampleSafety,
                                            G4double sampleStep) const
{
  if ( fVerbose == 1 )
  {
    G4cout << "Daughter "
           << std::setw(15) << sampleSafety << " ";
    if (sampleStep)
    {
      G4cout << std::setw(15) << sampleStep << " ";
    }
    else
    {
      G4cout << std::setw(15) << "N/C" << " ";
    }
    G4cout << samplePoint   << " - "
           << sampleSolid->GetEntityType() << ": " << sampleSolid->GetName()
           << G4endl;
  }
}
