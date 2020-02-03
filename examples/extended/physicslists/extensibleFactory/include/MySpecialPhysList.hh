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
//
//---------------------------------------------------------------------------
//
// ClassName:  MySpecialPhysList
//
// Author: 2017-11-02 R. Hatcher
//   Example "alternative" physics list by typedef'ing QBBC
//
//----------------------------------------------------------------------------
//
#ifndef TMySpecialPhysList_h
#define TMySpecialPhysList_h 1

#include "MySpecialPhysList.icc"

// Users would define their own physics list in
//   MySpecialPhysList.icc and MySpecialPhysList.hh
//
// The only requirement for registering with the extensible factory
// is that the physics list constructor must accept a single G4int argument
// which is the verbosity, i.e. :-
//
//  class MySpecialPhysicsList : public G4ModularPhysList
//  {
//     public:
//       MySpecialPhysicsList(G4int ver = 1 [, any defaulted args ] );
//       virtual ~MySpecialPhysicsList();
//       virtual void SetCuts();
//     ....

#include "QBBC.hh"
typedef   QBBC MySpecialPhysList;

#endif

