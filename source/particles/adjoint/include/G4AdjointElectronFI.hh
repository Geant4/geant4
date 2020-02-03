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
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//    1st March 2007 creation by L. Desorgher based on a modification of G4Electron 	   
//    06  Nov.  2008 modified for Geant4-09-02  by Hisaya Kurashige 
//
//-------------------------------------------------------------
//	Documentation:
//   Adjoint particles are used in Reverse/Adjoint Monte Carlo simulations. New adjoint 
//   processes act on adjoint particles when they are  tracked backward in the geometry. 
//   The use of adjoint particles instead of "normal" particles 
//   during a reverse simulation is based on an idea of M. Asai   
//
#ifndef G4AdjointElectronFI_h
#define G4AdjointElectronFI_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                         ADJOINT ELECTRON                       ###
// ######################################################################

class G4AdjointElectronFI : public G4ParticleDefinition
{
 private:
   static G4AdjointElectronFI* theInstance;
   G4AdjointElectronFI(){}
   ~G4AdjointElectronFI(){}

 public:
   static G4AdjointElectronFI* Definition();
   static G4AdjointElectronFI* AdjointElectronFIDefinition();
   static G4AdjointElectronFI* AdjointElectronFI();
};

#endif











