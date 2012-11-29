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
// $Id$
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Quasmon initiated Candidates used by CHIPS model
// -----------------------------------------------------------------
// Short description: A candidate for hadronization. The candidates
// (hadrons or nuclear fragments) are competative, each quark of a
// Quasmon select which candidate to use for hadronization
// ------------------------------------------------------------------

#ifndef G4QCandidateVector_h
#define G4QCandidateVector_h 1

#include "G4QCandidate.hh"
#include <vector>

typedef std::vector<G4QCandidate *> G4QCandidateVector;
struct DeleteQCandidate{ void operator()(G4QCandidate * aQ){delete aQ;}};

#endif
