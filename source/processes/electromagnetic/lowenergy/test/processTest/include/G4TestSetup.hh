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
// $Id: G4TestSetup.hh,v 1.8 2006-06-29 19:48:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test DoIt method of physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#ifndef G4TESTSETUP_HH
#define G4TESTSETUP_HH 1

#include "globals.hh"

class G4Material;
class G4ParticleDefinition;
class G4PVPlacement;
class G4Track;
class G4Step;

class G4TestSetup {
 
public:

  G4TestSetup(const G4Material* aMaterial,
	      const G4ParticleDefinition* def,
	      G4double minEnergy, G4double  maxEnergy);

  virtual ~G4TestSetup();
 
  void makeGeometry();
  const G4Track* makeTrack();
  const G4Step* makeStep();

private:
  
  // Hide copy constructor and assignment operator
  G4TestSetup(const G4TestSetup&);
  G4TestSetup& operator=(const G4TestSetup& right);

  G4ParticleDefinition* part;
  G4Material* material;
  G4double eMin;
  G4double eMax;
  G4Track* track;
  G4Step* step;
  G4PVPlacement* physicalFrame;

};
 
#endif
