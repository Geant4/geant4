// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRFeatures.hh,v 1.2 1999-01-09 16:11:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#if !defined G4_FR_FEATURES_HH
#define      G4_FR_FEATURES_HH


const char FR_DAWN_FEATURES[] = "High quality technical renderer.\
\n    Features:      exact hidden line, hidden surface algorithms.\
\n                   high (unlimited) resolution.\
\n                   renders to PostScript for viewing and/or hardcopy.\
\n                   remote rendering.\
\n                   off-line rendering.\
\n                   graphical user interface.\
\n                   connection via TCP/IP socket or named pipe Fukui Renderer DAWN.\
\n    Disadvantages: compute intensive, takes time (use a fast graphics\
\n                   system, such as OpenGL, to select view, then copy\
\n                   to this renderer - /vis~/copy/view, /vis~/set/view).";


const char FR_DAWNFILE_FEATURES[] ="High quality technical renderer.\
\n    Features:      exact hidden line, hidden surface algorithms.\
\n                   high (unlimited) resolution.\
\n                   renders to PostScript for viewing and/or hardcopy.\
\n                   remote rendering.\
\n                   off-line rendering.\
\n                   graphical user interface.\
\n                   connection via g4.prim file to Fukui Renderer DAWN,\
\n                   DAVID (DAwn's Visual Intersection Debugger, etc.\
\n    Disadvantages: compute intensive, takes time (use a fast graphics\
\n                   system, such as OpenGL, to select view, then copy\
\n                   to this renderer - /vis~/copy/view, /vis~/set/view).";

#endif
