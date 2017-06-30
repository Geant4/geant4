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
// $Id: G4RayTracerSceneHandler.hh 104015 2017-05-08 07:28:08Z gcosmo $

// John Allison  17th March 2000

#ifndef G4RAYTRACERSCENEHANDLER_HH
#define G4RAYTRACERSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

#include "G4ModelingParameters.hh"
#include <map>

class G4RayTracerSceneHandler: public G4VSceneHandler
{
public:

  G4RayTracerSceneHandler(G4VGraphicsSystem& system,
			  const G4String& name = "");
  virtual ~G4RayTracerSceneHandler();

  void BuildVisAttsMap (const G4VSolid& solid);
  // BuildVisAttsMap is the appropriate name for what it does in this class
  // but it is actually an override of G4VSceneHandler::RequestPrimitives.
  void RequestPrimitives (const G4VSolid& solid)
  {BuildVisAttsMap(solid);}

  // Required pure virtual functions, not used in Ray Tracer.
  void AddPrimitive(const G4Polyline&){}
  void AddPrimitive(const G4Text&){}
  void AddPrimitive(const G4Circle&){}
  void AddPrimitive(const G4Square&){}
  void AddPrimitive(const G4Polyhedron&){}
  void AddPrimitive(const G4Polymarker&){}
  void AddPrimitive(const G4Scale&){}

  void ClearStore ();

  struct PathLessThan
  {G4bool operator()
    (const G4ModelingParameters::PVPointerCopyNoPath&,
     const G4ModelingParameters::PVPointerCopyNoPath&) const;
  };

  const std::map
  <G4ModelingParameters::PVPointerCopyNoPath,G4VisAttributes,PathLessThan>&
  GetSceneVisAttsMap() const
  {return fSceneVisAttsMap;}

private:
  
  static G4int fSceneIdCount;  // Counter for RayTracer scene handlers.

  std::map
  <G4ModelingParameters::PVPointerCopyNoPath,G4VisAttributes,PathLessThan>
  fSceneVisAttsMap;
};

#endif
