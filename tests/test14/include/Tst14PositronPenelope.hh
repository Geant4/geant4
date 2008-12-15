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
// $Id: Tst14PositronPenelope.hh,v 1.1 2008-12-15 10:23:20 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola (pandola@lngs.infn.it)
//
// History:
// -----------
// 15 Dec 2008 Luciano Pandola     Created
//
// -------------------------------------------------------------------

// Class description:
// System test for e/gamma, positron processes a' la Penelope for PhysicsList

// -------------------------------------------------------------------

#ifndef TST14POSITRONPENELOPE_HH
#define TST14POSITRONPENELOPE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst14PositronPenelope : public G4VPhysicsConstructor {

public: 

  Tst14PositronPenelope(const G4String& name = "positron-penelope");
  
  virtual ~Tst14PositronPenelope();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif








