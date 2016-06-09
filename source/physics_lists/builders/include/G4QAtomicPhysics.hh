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
// $Id: G4QAtomicPhysics.hh,v 1.2 2010-06-03 14:37:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QAtomicPhysics
//
// Author:      M. Kosov 20.11.2009 (similar to G4EmStandardPhysics)
//
// Modified:
//
//---------------------------------------------------------------------------
//
// This class provides construction of CHIPS-modified Electromagnetic physics
//
//---------------------------------------------------------------------------

#ifndef G4QAtomicPhysics_h
#define G4QAtomicPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class G4QAtomicPhysics : public G4VPhysicsConstructor
{
public:
  G4QAtomicPhysics(G4int ver = 0);
  G4QAtomicPhysics(G4int ver, const G4String& name);
  virtual ~G4QAtomicPhysics();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
};

#endif






