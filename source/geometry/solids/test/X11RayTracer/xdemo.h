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
unsigned char        *buffer;

GC                   gc;
XGCValues            gcVal;
unsigned long        gcMask;

/* */

void init_x ();
int  event_x ();
void close_x ();

#endif
