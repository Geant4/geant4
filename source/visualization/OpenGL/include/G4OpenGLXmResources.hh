// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmResources.hh,v 1.1 1999-01-07 16:14:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Default resources file for GEANT4 OpenGL Motif windows.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMRESOURCES_HH
#define G4OPENGLXMRESOURCES_HH

static String fallbackResources[] = {
  "*glxarea*width: 500", "*glxarea*height: 500",
  "*frame*x: 10", "*frame*y: 10",
  "*frame*topOffset: 10", "*frame*bottomOffset: 10",
  "*frame*rightOffset: 10", "*frame*leftOffset: 10",
  "*frame*shadowType: SHADOW_IN", "*useColorObj: False", 
  NULL
  };

#endif

#endif
