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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP_INCLXX
//
// Author: 2011 P. Kaitaniemi
//
// Modified:
// D. Mancusi 23.03.2012 Extended INCL++ to incident heavy ions up to 16O
//
//----------------------------------------------------------------------------
//
#ifndef TQGSP_INCLXX_h
#define TQGSP_INCLXX_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

/**
 * <h1>Physics list QGSP_INCLXX</h1>
 *
 * <h2>Use case</h2>
 * This list is mainly intended for use with energies less than 3
 * GeV. This is useful for e.g. spallation studies and Accelerator
 * Driven Systems (ADS) applications.
 *
 * <h2>Usage</h2>
 * The physics list can be activated in a simulation application by
 * giving it as part of the user initialization to the run manager:
 * @code
 * G4RunManager *runManager = new G4RunManager;
 * G4VUserPhysicsList *physics = new QGSP_INCLXX;
 * runManager->SetUserInitialization(physics);
 * @endcode
 *
 * <h2>Hadronic models</h2>
 * The list uses INCL++ intra-nuclear cascade model in the energy range 0 - 3
 * GeV. Between 3 - 15 GeV Bertini cascade is used and above 15 GeV the high
 * energy QGSP model.
 *
 * @see HadronPhysicsQGSP_INCLXX
 * @see G4INCLXXProtonBuilder
 * @see G4INCLXXNeutronBuilder
 * @see G4INCLXXPiKBuilder
 * @see G4IonINCLXXPhysics
 */
template<class T>
class TQGSP_INCLXX: public T
{
public:
  TQGSP_INCLXX(G4int ver = 1);
  virtual ~TQGSP_INCLXX();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};

#include "QGSP_INCLXX.icc"
typedef TQGSP_INCLXX<G4VModularPhysicsList> QGSP_INCLXX;

// 2011 by P. Kaitaniemi

#endif



