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
// $Id: G4OpenGLXmResources.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
