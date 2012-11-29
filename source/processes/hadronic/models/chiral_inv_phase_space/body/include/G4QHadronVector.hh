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
// Type defenition for a Vector of Hadrons - output of CHIPS model
// ---------------------------------------------------------------
// Short description: In CHIPS all particles are G4QHadrons, while they
// can be leptons, gammas or nuclei. The G4QPDGCode makes the difference.
// In addition the 4-momentum is a basic value, so the mass can be
// different from the GS mass (e.g. for the virtual gamma). This class
// is made for the output list of hadrons.
// -------------------------------------------------------------------

#ifndef G4QHadronVector_h
#define G4QHadronVector_h 1
//
// $Id$
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4QParton ----------------
//             by Mikhail Kossov, Oct 1999.
// class for QHadronVector (string) used by CHIPS Models
// ------------------------------------------------------------

#include "G4QHadron.hh"
#include <vector>

typedef std::vector<G4QHadron *> G4QHadronVector;
struct DeleteQHadron { void operator()(G4QHadron* aQH){delete aQH;}};

#endif
