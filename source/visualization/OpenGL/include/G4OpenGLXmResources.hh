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
// $Id: G4OpenGLXmResources.hh,v 1.4 2001-07-11 10:08:51 gunter Exp $
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
