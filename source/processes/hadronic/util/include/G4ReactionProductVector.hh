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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ReactionProductVector.hh,v 1.6 2001-10-04 20:00:44 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	History: first implementation, alternative to G4FastVector
//               less fast, but it has variable array length and checks boundaries
//	26th September, Chr. Voelcker
// ------------------------------------------------------------

#ifndef G4ReactionProductVector_h
#define G4ReactionProductVector_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4ReactionProduct.hh"
#include "g4std/vector"

// #ifdef STL
// //in future use STL vector as container of reaction products ...
// typedef Vector<G4ReactionProduct> G4ReactionProductVector;
// #elseifdef RWT

typedef G4std::vector<G4ReactionProduct *> G4ReactionProductVector;
struct DeleteReactionProduct{ void operator()(G4ReactionProduct * aR){delete aR;} };

// #endif

#endif
