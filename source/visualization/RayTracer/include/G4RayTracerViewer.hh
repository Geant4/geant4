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
// $Id: G4RayTracerViewer.hh 103626 2017-04-19 13:29:18Z gcosmo $

// John Allison  17th March 2000

#ifndef G4RAYTRACERVIEWER_HH
#define G4RAYTRACERVIEWER_HH

#include "G4VViewer.hh"

class G4TheRayTracer;

class G4RayTracerViewer: public G4VViewer {
public:
  G4RayTracerViewer(G4VSceneHandler&,
		    const G4String& name,
		    G4TheRayTracer* = 0);
  virtual ~G4RayTracerViewer();
  void Initialise();
  void SetView();
  void ClearView();
  void DrawView();
  G4TheRayTracer* GetTracer() {return theTracer;}
protected:
  G4int fFileCount;
  G4TheRayTracer* theTracer;
};

#endif
