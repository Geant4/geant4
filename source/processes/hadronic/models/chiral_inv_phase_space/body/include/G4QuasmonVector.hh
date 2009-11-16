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
// $Id: G4QuasmonVector.hh,v 1.22 2009-11-16 18:15:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QuasmonVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type definition for a Vector of Quasmons - output of CHIPS model
// ----------------------------------------------------------------
// Short description: If partons from the G4QPartonPair are close in
// rapidity, they create Quasmons, but if they are far in the rapidity
// space, they can not interact directly. Say the bottom parton (quark)
// has rapidity 0, and the top parton (antiquark) has rapidity 8, then
// the top quark splits in two by radiating gluon, and each part has
// rapidity 4, then the gluon splits in quark-antiquark pair (rapidity
// 2 each), and then the quark gadiates anothe gluon and reachs rapidity
// 1. Now it can interact with the bottom antiquark, creating a Quasmon
// or a hadron. The intermediate partons is the string ladder.
// ---------------------------------------------------------------------

#ifndef G4QuasmonVector_h
#define G4QuasmonVector_h 1

#include "G4Quasmon.hh"
#include <vector>

typedef std::vector<G4Quasmon *> G4QuasmonVector;
struct DeleteQuasmon{ void operator()(G4Quasmon *aN){delete aN;} };

#endif
