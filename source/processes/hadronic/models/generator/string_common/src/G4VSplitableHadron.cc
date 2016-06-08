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
// $Id: G4VSplitableHadron.cc,v 1.6 2001/10/05 16:17:16 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
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
      :  theDefinition(NULL), theCollisionCount(0), isSplit(false)
{
}

G4VSplitableHadron::G4VSplitableHadron(const G4ReactionProduct & aPrimary)
      :  theCollisionCount(0), isSplit(false)
{
	theDefinition=aPrimary.GetDefinition();
	the4Momentum.setVect(aPrimary.GetMomentum());
	the4Momentum.setE(aPrimary.GetTotalEnergy());
}

G4VSplitableHadron::G4VSplitableHadron(const G4Nucleon & aNucleon)
{
	theCollisionCount=0;
  isSplit = false;
	theDefinition=aNucleon.GetParticleType();
	the4Momentum=aNucleon.GetMomentum();
	thePosition=aNucleon.GetPosition();
}

G4VSplitableHadron::G4VSplitableHadron(const G4VKineticNucleon * aNucleon)
{
	theCollisionCount=0;
  isSplit = false;
	theDefinition=aNucleon->GetDefinition();
	the4Momentum=aNucleon->Get4Momentum();
	thePosition=aNucleon->GetPosition();
}

G4VSplitableHadron::G4VSplitableHadron(const G4VSplitableHadron &right)
{
	theCollisionCount=0;
  isSplit = false;
	theDefinition= right.GetDefinition();
	the4Momentum= right.Get4Momentum();
	thePosition=  right.GetPosition();
}


G4VSplitableHadron::~G4VSplitableHadron()
{
}


const G4VSplitableHadron & G4VSplitableHadron::operator=(const G4VSplitableHadron &right)
{
  G4Exception("G4VSplitableHadron::operator= meant to not be accessable");
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
