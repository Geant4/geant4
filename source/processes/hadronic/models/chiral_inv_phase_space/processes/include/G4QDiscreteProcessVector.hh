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
//      ---------------- G4QDiscreteProcessVector ----------------
//             by Mikhail Kossov, Aug 2007.
// Type defenition for Vectors of DiscreteProcesses (G4VDiscreteProcess)
// ---------------------------------------------------------------------
// Short description: universal mixer of processes (NOT models as in GHAD!)
// depending on the application energy region G4QDiscProcessMixer is using
// this vector of the Geant4 processes.
// ------------------------------------------------------------------------

#ifndef G4QDiscreteProcessVector_h
#define G4QDiscreteProcessVector_h 1

#include "G4VDiscreteProcess.hh"
#include <vector>

typedef std::vector<std::pair<G4VDiscreteProcess*, G4double>*> G4QDiscreteProcessVector;
struct DeleteDiscreteProcess { void operator()(std::pair<G4VDiscreteProcess*,G4double>* DP)
                               {delete DP->first; delete DP;} };

#endif
