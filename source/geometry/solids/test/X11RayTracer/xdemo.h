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
