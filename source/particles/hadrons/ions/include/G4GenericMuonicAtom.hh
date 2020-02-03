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
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, K.L. Genser 2 December 2016
//      based on G4GenericIon
//
// ****************************************************************
// This class is used only by G4IonTable and not for tracking
// One should register processes for G4MuonicAtoms with this class
//
// ----------------------------------------------------------------------

#ifndef G4GenericMuonicAtom_h
#define G4GenericMuonicAtom_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonicAtom.hh"

// ######################################################################
// ###                          GenericMuonicAtom                     ###
// ######################################################################

class G4GenericMuonicAtom : public G4MuonicAtom
{
 private:
   static G4GenericMuonicAtom* theInstance;
   G4GenericMuonicAtom(){}
   ~G4GenericMuonicAtom(){}

 public:
   static G4GenericMuonicAtom* Definition();
   static G4GenericMuonicAtom* GenericMuonicAtomDefinition();
   static G4GenericMuonicAtom* GenericMuonicAtom();
};

#endif
