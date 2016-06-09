//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RayTracer.hh,v 1.10 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $
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

class G4RayTracer : public G4VGraphicsSystem
{
  public: // with description
    G4RayTracer();
    ~G4RayTracer();
    G4VSceneHandler* CreateSceneHandler (const G4String& );
    G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& );
};

#endif
