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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NeutronHPCaptureFS.hh,v 1.8 2002-12-12 19:18:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPCaptureFS_h
#define G4NeutronHPCaptureFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4ReactionProductVector.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPPhotonDist.hh"

class G4NeutronHPCaptureFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPCaptureFS()
  {
    hasXsec = false; 
    targetMass = 0;
  }
  
  ~G4NeutronHPCaptureFS()
  {
  }
  
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPCaptureFS * theNew = new G4NeutronHPCaptureFS;
   return theNew;
  }
  
  private:
  
  G4double targetMass;
  
  G4NeutronHPPhotonDist theFinalStatePhotons;
  
  G4NeutronHPNames theNames;
  
  G4double theCurrentA;
  G4double theCurrentZ;
};
#endif
