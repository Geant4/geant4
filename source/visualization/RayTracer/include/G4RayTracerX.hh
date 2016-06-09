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
// $Id: G4RayTracerX.hh,v 1.2 2005/07/20 20:39:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// John Allison  11th September 2003

#if defined (G4VIS_BUILD_RAYTRACERX_DRIVER) || (G4VIS_USE_RAYTRACERX)

#ifndef G4RAYTRACERX_HH
#define G4RAYTRACERX_HH

// class description:
//
// G4RayTracerX
// X window version of RayTracer (also produces jpeg file - see class
// description for G4RayTracer).

#include "G4RayTracer.hh"

#include <X11/Xlib.h>
#include <X11/Xutil.h>

class G4RayTracerX : public G4RayTracer
{
public:
  G4RayTracerX();
  G4VSceneHandler* CreateSceneHandler (const G4String& );
  G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& );

  // X Window variables...
  Display* display;
  Window win;
  GC gc;
  XStandardColormap *scmap;
};

#endif

#endif
