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
// $Id: G4NeutronHPLCFissionFS.hh,v 1.6 2002-12-12 19:18:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPLCFissionFS_h
#define G4NeutronHPLCFissionFS_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFissionBaseFS.hh"

class G4NeutronHPLCFissionFS : public G4NeutronHPFissionBaseFS
{
  public:
  
  G4NeutronHPLCFissionFS(){ hasXsec = false; }
  ~G4NeutronHPLCFissionFS(){}
  void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType);
  G4DynamicParticleVector * ApplyYourself(G4int NNeutrons);
  G4NeutronHPFinalState * New() 
  {
   G4NeutronHPLCFissionFS * theNew = new G4NeutronHPLCFissionFS;
   return theNew;
  }
  
  private:
  G4ParticleChange * ApplyYourself(const G4Track & theTrack) { return NULL; }
    
};
#endif
