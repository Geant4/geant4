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
// $Id: G4QPartonPairVector.hh,v 1.2 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, October 2006.
// Type defenition for a Vector of PartonPairs of CHIPS model
// ---------------------------------------------------------------
// Short description: Each Quasmon String has a pair of partons
// (quark/diquark-partons) on its ends. During the hadronization
// procedure the rapidity gap between partons shrinks, but the
// parton pair still exists, while it is converted to the final
// meson (quaek-antiquark) or baryon (quark-diquark).
// --------------------------------------------------------------

#ifndef G4QPartonPairVector_h
#define G4QPartonPairVector_h 1

#include "G4QPartonPair.hh"
#include <vector>

typedef std::vector<G4QPartonPair *> G4QPartonPairVector;
struct DeleteQPartonPair { void operator()(G4QPartonPair* aQPP){delete aQPP;}};

#endif
