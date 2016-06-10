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
// $Id: G4UnknownParticle.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//  first implementation : M.Asai Jul 07, 2004
// ----------------------------------------------------------------
//  New impelemenataion as an utility class  H.Kurashige, 14 July 2004
// ----------------------------------------------------------------


#ifndef G4UnknownParticle_h
#define G4UnknownParticle_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                         UNKNOWN                                ###
// ######################################################################

class G4UnknownParticle : public G4ParticleDefinition
{
 private:
   static G4UnknownParticle* theInstance;

 private:
  G4UnknownParticle(){}

 public:
   ~G4UnknownParticle(){}
 
   static G4UnknownParticle* Definition();
   static G4UnknownParticle* UnknownParticleDefinition();
   static G4UnknownParticle* UnknownParticle();
};

#endif











