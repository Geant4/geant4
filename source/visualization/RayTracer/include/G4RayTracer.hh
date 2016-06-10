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
// $Id: G4RayTracer.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
//


#ifndef G4RayTracer_H
#define G4RayTracer_H 1

// class description:
//
// G4RayTracer
// This is a factory for G4RayTracerSceneHandler and G4RayTracerViewer
// objects.  The viewer uses the original G4RayTracer, which has now
// been renamed G4TheRayTracer.

#include "G4VGraphicsSystem.hh"

class G4TheRayTracer;

class G4RayTracer : public G4VGraphicsSystem
{
public: // with description
  G4RayTracer();
  ~G4RayTracer();
  G4VSceneHandler* CreateSceneHandler (const G4String& );
  G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& );
private:
  G4TheRayTracer* theRayTracer;
};

#endif
