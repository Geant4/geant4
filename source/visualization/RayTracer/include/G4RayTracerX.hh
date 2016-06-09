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
// $Id: G4RayTracerX.hh,v 1.3 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $
//
//
// John Allison  11th September 2003

#if defined (G4VIS_BUILD_RAYTRACERX_DRIVER) || (G4VIS_USE_RAYTRACERX)

#ifndef G4RAYTRACERX_HH
#define G4RAYTRACERX_HH

// class description:
//
// G4RayTracerX
// X window version of RayTracer - opens X window for viewing as well
// as producing jpeg file just like G4RayTracer.

#include "G4VGraphicsSystem.hh"

class G4RayTracerX : public G4VGraphicsSystem
{
  public: // with description
    G4RayTracerX();
    ~G4RayTracerX();
    G4VSceneHandler* CreateSceneHandler (const G4String& );
    G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& );
};

#endif
#endif
