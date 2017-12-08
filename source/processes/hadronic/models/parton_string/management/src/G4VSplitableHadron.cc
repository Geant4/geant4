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
// $Id: G4VSplitableHadron.cc 107318 2017-11-08 16:27:32Z gcosmo $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4VSplitableHadron----------------
//             by Gunter Folger, June 1998.
//       class storing an interacting particle. Used by Parton String Models.
// ------------------------------------------------------------

#include "G4VSplitableHadron.hh"
#include "G4Nucleon.hh"
#include "G4VKineticNucleon.hh"

G4VSplitableHadron::G4VSplitableHadron()
: theDefinition(NULL), TimeOfCreation(0.), theCollisionCount(0), curStatus(0), isSplit(false)
{
}

G4VSplitableHadron::G4VSplitableHadron(const G4ReactionProduct & aPrimary)
: TimeOfCreation(0.), theCollisionCount(0), curStatus(0), isSplit(false)
{
  theDefinition=aPrimary.GetDefinition();
  the4Momentum.setVect(aPrimary.GetMomentum());
  the4Momentum.setE(aPrimary.GetTotalEnergy());
}

G4VSplitableHadron::G4VSplitableHadron(const G4Nucleon & aNucleon)
{
  TimeOfCreation   = 0.;                                          
  theCollisionCount= 0;
  isSplit          = false;
  theDefinition    = aNucleon.GetParticleType();
  the4Momentum     = aNucleon.GetMomentum();
  thePosition      = aNucleon.GetPosition();
  curStatus        = 0;                                           
}

G4VSplitableHadron::G4VSplitableHadron(const G4VKineticNucleon * aNucleon)
{
  TimeOfCreation   = 0.;   
  theCollisionCount= 0;
  isSplit          = false;
  theDefinition    = aNucleon->GetDefinition();
  the4Momentum     = aNucleon->Get4Momentum();
  thePosition      = aNucleon->GetPosition();
  curStatus        = 0;                                        
}

G4VSplitableHadron::G4VSplitableHadron(const G4VSplitableHadron &right)
{
  TimeOfCreation   = 0.;
  theCollisionCount= 0;
  isSplit          = false;
  theDefinition    = right.GetDefinition();
  the4Momentum     = right.Get4Momentum();
  thePosition      = right.GetPosition();
  curStatus        = 0;                                        
}


G4VSplitableHadron::~G4VSplitableHadron()
{
}


const G4VSplitableHadron & G4VSplitableHadron::operator=(const G4VSplitableHadron &)
{
  throw G4HadronicException(__FILE__, __LINE__, 
                            "G4VSplitableHadron::operator= meant to not be accessable");
  return *this;
}


int G4VSplitableHadron::operator==(const G4VSplitableHadron &right) const
{
  return this==&right;
}

int G4VSplitableHadron::operator!=(const G4VSplitableHadron &right) const
{
  return this!=&right;
}


void G4VSplitableHadron::SplitUp()
{
}

