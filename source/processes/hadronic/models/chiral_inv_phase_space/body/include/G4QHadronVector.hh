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
//
// $Id: G4QHadronVector.hh,v 1.11 2002-12-12 19:14:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for a Vector of Hadrons - output of CHIPS model
// ---------------------------------------------------------------

#ifndef G4QHadronVector_h
#define G4QHadronVector_h 1

#include "G4QHadron.hh"
#include "g4std/vector"

typedef G4std::vector<G4QHadron *> G4QHadronVector;
struct DeleteQHadron { void operator()(G4QHadron * aQ){delete aQ;}};

#endif
