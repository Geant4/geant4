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
// $Id: G4FRFeatures.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
