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
// $Id: G4NeutronHPFissionFS.hh,v 1.5 2001-07-26 09:28:02 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFissionFS_h
#define G4NeutronHPFissionFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"

#include "G4NeutronHPFCFissionFS.hh"
#include "G4NeutronHPSCFissionFS.hh"
#include "G4NeutronHPTCFissionFS.hh"
#include "G4NeutronHPLCFissionFS.hh"
#include "G4NeutronHPFSFissionFS.hh"

class G4NeutronHPFissionFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFissionFS(){ hasXsec = false; }
  ~G4NeutronHPFissionFS(){}
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4ParticleChange * ApplyYourself(const G4Track & theTrack);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPFissionFS * theNew = new G4NeutronHPFissionFS;
   return theNew;
  }
        
  private:
  
  G4NeutronHPFSFissionFS theFS;
  G4NeutronHPFCFissionFS theFC;
  G4NeutronHPSCFissionFS theSC;
  G4NeutronHPTCFissionFS theTC;
  G4NeutronHPLCFissionFS theLC;
    
};
#endif
