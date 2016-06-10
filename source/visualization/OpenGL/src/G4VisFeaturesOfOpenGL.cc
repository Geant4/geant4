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
// $Id: G4VisFeaturesOfOpenGL.cc 66264 2012-12-14 10:17:44Z allison $
//
// 

#include "G4VisFeaturesOfOpenGL.hh"

G4String G4VisFeaturesOfOpenGLIX () {
  return
    "    Dumb single buffered X Window, No Graphics Database."
    "\n    Advantages:    does not gobble server memory."
    "\n                   good for drawing steps and hits."
    "\n    Disadvantages: needs G4 kernel for re-Draw."
    "\n                   cannot take advantage of graphics accelerators.";
}

G4String G4VisFeaturesOfOpenGLSX () {
  return
    "    Dumb double buffered X Window with Graphics Database."
    "\n    Advantages:    uses display lists as graphics database."
    "\n                   fastest possible redraw, e.g., on simple change"
    "\n                     of viewpoint."
    "\n                   uses client-server model for remote viewing"
    "\n                     (but only if you have a full client-server"
    "\n                     implementation of OpenGL, i.e., not Mesa)."
    "\n    Disadvantages: not advised for viewing large numbers of steps"
    "\n                     and/or hits, because it gobbles memory for"
    " database.";
}

G4String G4VisFeaturesOfOpenGLIXm () {
  return
    "    Smart single buffered X Window, No Graphics Database."
    "\n    Advantages:    resizeable, and has Motif-based view-control panel."
    "\n                   does not gobble server memory."
    "\n                   good for drawing steps and hits."
    "\n    Disadvantages: currently locks out GEANT4 commands, until \"exit\"."
    "\n                   needs G4 kernel for re-Draw."
    "\n                   cannot take advantage of graphics accelerators.";
}

G4String G4VisFeaturesOfOpenGLSXm () {
  return
    "    Smart double buffered X Window with Graphics Database."
    "\n    Advantages:    resizeable, and has Motif-based view-control panel."
    "\n                   uses display lists as graphics database."
    "\n                   fastest possible redraw, e.g., on simple change"
    "\n                     of viewpoint."
    "\n                   uses client-server model for remote viewing"
    "\n                     (but only if you have a full client-server"
    "\n                     implementation of OpenGL, i.e., not Mesa)."
    "\n    Disadvantages: currently locks out GEANT4 commands, until \"exit\"."
    "\n                   not advised for viewing large numbers of steps"
    "\n                     and/or hits, because it gobbles memory for"
    " database.";
}

G4String G4VisFeaturesOfOpenGLIWin32 () {
  return
    "\n It runs on WindowsNT ";
}

G4String G4VisFeaturesOfOpenGLSWin32 () {
  return
    "\n It runs on WindowsNT ";
}

G4String G4VisFeaturesOfOpenGLIQt () {
  return
    "\n It runs everywhere";
}

G4String G4VisFeaturesOfOpenGLSQt () {
  return
    "\n It runs everywhere ";
}

G4String G4VisFeaturesOfOpenGLIWt () {
  return
    "\n It runs everywhere";
}
