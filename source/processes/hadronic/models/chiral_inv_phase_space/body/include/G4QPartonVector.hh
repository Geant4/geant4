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
//             by Mikhail Kossov, Oct 2006.
// Type defenition for a Vector of Partons - output of CHIPS model
// ---------------------------------------------------------------
// Short description: The Quark-Gluon String consists of the partons, which
// are quarks and some times gluons.
// ------------------------------------------------------------------------

#ifndef G4QPartonVector_h
#define G4QPartonVector_h 1
//
// $Id$
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QParton ----------------
//             by Mikhail Kossov, Oct 2006.
// class for PartonVector (string) used by Parton String Models
// ------------------------------------------------------------

#include "G4QParton.hh"
#include <vector>

typedef std::vector<G4QParton *> G4QPartonVector;
struct DeleteQParton { void operator()(G4QParton* aQP){delete aQP;}};

#endif
