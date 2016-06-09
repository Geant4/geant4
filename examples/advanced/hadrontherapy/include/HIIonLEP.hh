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
// $Id: HadrontherapyProtonPrecompound.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

// ==============================
// PHYSICS PROCESSES:
// ==============================
//  Hadronic inelastic processes for: 
//      deuterons, tritons, alpha particles
//     
// 
// ==============================
// COMMENTS:
// ==============================
//  The considered process is G4XXXInelasticProcess (XXX to be replaced
//  by the particular ion name), and the applied model is G4LEXXXInelastic
//  Considered cross sections are G4TripathiCrossSection and 
//  G4IonsShenCrossSection
//

#ifndef HIIONLEP_H
#define HIIONLEP_H 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh" 

class HIIonLEP: public G4VPhysicsConstructor 
{
 public:
   HIIonLEP(const G4String& name = "Ion-LEP");
   virtual ~HIIonLEP();

 protected:
   void ConstructParticle(){};
   void ConstructProcess();
};
#endif 








