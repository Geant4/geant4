//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPCapture.hh,v 1.6 2001-07-26 09:27:54 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron capture below 20 MeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#ifndef G4NeutronHPCapture_h
#define G4NeutronHPCapture_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"

class G4NeutronHPCapture : public G4HadronicInteraction
{
  public: 
  
  G4NeutronHPCapture();

  ~G4NeutronHPCapture();

  G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
  
  G4double * xSec;
  G4NeutronHPChannel * theCapture;
  G4String dirName;
  G4int numEle;
  G4int it;
  
  G4ParticleChange theResult;
};

#endif
