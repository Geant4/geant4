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
// $Id$
//
#ifndef G4NeutronHPFissionFS_h
#define G4NeutronHPFissionFS_h 1

#include "globals.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"

#include "G4NeutronHPFCFissionFS.hh"
#include "G4NeutronHPSCFissionFS.hh"
#include "G4NeutronHPTCFissionFS.hh"
#include "G4NeutronHPLCFissionFS.hh"
#include "G4NeutronHPFSFissionFS.hh"

#include "G4NeutronHPFFFissionFS.hh"

class G4NeutronHPFissionFS : public G4NeutronHPFinalState
{
  public:
  
  G4NeutronHPFissionFS(){ hasXsec = false; produceFissionFragments = false; }
  ~G4NeutronHPFissionFS(){}
  void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType);
  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack);
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
    
  G4NeutronHPFFFissionFS theFF;
  G4bool produceFissionFragments;
};
#endif
