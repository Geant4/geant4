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
// $Id: G4NeutronHPFission.hh,v 1.10 2006-06-29 20:47:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
  
#ifndef G4NeutronHPFission_h
#define G4NeutronHPFission_h 1

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron induced fission below 20 MeV; 
// Note that this model (by intent of avoiding the possibility of heating studies) does
// not provide the nuclear fragments.
//
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

#include "G4NeutronHPFissionFS.hh"

class G4NeutronHPFission : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPFission();

  ~G4NeutronHPFission();
  
  G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);

  private:
  
  G4NeutronHPFissionFS theFS;
  
  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theFission;
  G4String dirName;
  G4int numEle;
  // static G4String theNames[3];
};

#endif
