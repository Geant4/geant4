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
// $Id: G4VisFeaturesOfOpenGL.cc,v 1.4 2001-07-11 10:09:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
