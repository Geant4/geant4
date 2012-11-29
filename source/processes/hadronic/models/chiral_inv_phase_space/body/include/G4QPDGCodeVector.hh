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
// Type defenition for Hadron definition in CHIPS model
// ---------------------------------------------------------------
// Short description: The PDG Code is made on the basis of the Quark
// Content (G4QuarkContent) of the hadronic state (including nuclear
// fragments). The PDG code of the ground state (e.g. pi, N, etc.) is
// calculated. It includes a complicated algortithm of the G.S. mass
// calculation for nuclear fragments (now it is synchronised with the
// G4 nuclear massess).
// -------------------------------------------------------------------

#ifndef G4QPDGCodeVector_h
#define G4QPDGCodeVector_h 1

#include "G4QPDGCode.hh"
#include <vector>

typedef std::vector<G4QPDGCode *> G4QPDGCodeVector;
struct DeleteQPDGCode
{
  void operator()(G4QPDGCode *aN)
  {
    //G4cout<<"G4QPDGCodeVector::DeleteQPDGCode: Before aN="<<aN<<G4endl; // TMP
    if(aN) delete aN; // void Destructor
    else G4cerr<<"***G4QPDGCodeVector::DeleteQPDGCode: aN="<<aN<<G4endl; // TMP
  }
};


#endif


