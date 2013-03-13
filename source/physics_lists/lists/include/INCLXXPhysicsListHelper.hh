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
// $Id: INCLXXPhysicsListHelper.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   INCLXXPhysicsListHelper
//
// Author: 04.03.2013 D. Mancusi
//   Created by modifying and generalising the QGSP_INCLXX source files
//
//----------------------------------------------------------------------------
//
#ifndef TINCLXXPhysicsListHelper_h
#define TINCLXXPhysicsListHelper_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"
#include "G4String.hh"

/**
 * <h1>Physics list helper INCLXXPhysicsListHelper</h1>
 *
 * <h2>Use case</h2>
 * This class is used to define the following INCLXX-based physics lists:
 * * QGSP_INCLXX
 * * QGSP_INCLXX_HP
 * * FTFP_INCLXX
 * * FTFP_INCLXX_HP
 * They are mainly intended for use with energies less than 3 GeV. This is
 * useful for e.g. spallation studies and Accelerator Driven Systems (ADS)
 * applications.
 *
 * <h2>Usage</h2>
 * The physics lists above can be activated in a simulation application by
 * giving it as part of the user initialization to the run manager, e.g.:
 * @code
 * G4RunManager *runManager = new G4RunManager;
 * G4VUserPhysicsList *physics = new QGSP_INCLXX;
 * runManager->SetUserInitialization(physics);
 * @endcode
 *
 * <h2>Hadronic models</h2>
 * The list uses INCL++ intra-nuclear cascade model in the energy range 0 - 3
 * GeV. Between 3 - 15 GeV Bertini cascade is used and above 15 GeV the high
 * energy QGSP model or FTFP model. The _HP variants use NeutronHP below 20
 * MeV.
 *
 * @see G4HadronPhysicsINCLXX
 * @see G4INCLXXProtonBuilder
 * @see G4INCLXXNeutronBuilder
 * @see G4INCLXXPiKBuilder
 * @see G4IonINCLXXPhysics
 */
template<class T, bool withNeutronHP, bool withFTFP>
class TINCLXXPhysicsListHelper: public T
{
public:
  TINCLXXPhysicsListHelper(G4int ver = 1);
  virtual ~TINCLXXPhysicsListHelper();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
  G4String name;
};

#include "INCLXXPhysicsListHelper.icc"
typedef TINCLXXPhysicsListHelper<G4VModularPhysicsList,false,false> QGSP_INCLXX;
typedef TINCLXXPhysicsListHelper<G4VModularPhysicsList,true, false> QGSP_INCLXX_HP;
typedef TINCLXXPhysicsListHelper<G4VModularPhysicsList,false,true> FTFP_INCLXX;
typedef TINCLXXPhysicsListHelper<G4VModularPhysicsList,true, true> FTFP_INCLXX_HP;

// 2013 by D. Mancusi

#endif



