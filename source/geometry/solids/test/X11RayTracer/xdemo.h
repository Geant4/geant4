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
#ifndef XDEMO_H
#define XDEMO_H

#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

/*
 *  Variables pour l'interface X
 */

#define X            0
#define Y            0
#define W            300
#define H            300
#define BORDER       0
#define NAME         "Emm G4 Visu EmmTorus"

Display              *dis;
int                  screen;
int                  depth;
int                  width;
int                  height;

Window               win;
Window               winRoot;
XSetWindowAttributes winAttr;
unsigned long        winMask;
XSizeHints           winHint;

XImage               *xim;
char                 *buffer;

GC                   gc;
XGCValues            gcVal;
unsigned long        gcMask;

/* */

void init_x ();
int  event_x ();
void close_x ();

#endif
