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
//---------------------------------------------------------------------------
//
// ClassName:      G4HadronicParameters
//
// Author:         2018 Alberto Ribon
//
// Description:    Singleton to keep global hadronic parameters.
//
//                 For the time being, at least, offers only "getters" but
//                 not "setters", i.e. a recompilation is needed to change
//                 the default parameters.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronicParameters_h
#define G4HadronicParameters_h 1

#include "globals.hh"


class G4HadronicParameters {
  public:

    static G4HadronicParameters* Instance();
    ~G4HadronicParameters();

    G4double GetMaxEnergy() const;
    // Upper limit for Geant4 hadronic physics, for any application.
    // Any hadronic model, physics list builder and constructor should use this method
    // instead of putting an arbitrary value in the code.
    // Any application which tries to use hadronic physics for an energy higher than this limit
    // will get a run-time crash, because no model is found.

  private:
    G4HadronicParameters();
    static G4HadronicParameters* sInstance;

    G4double fMaxEnergy;
};


inline G4double G4HadronicParameters::GetMaxEnergy() const { 
  return fMaxEnergy;
}

#endif

