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
// $Id: G4HEPEvtParticle.cc,v 1.3 2001-07-11 09:58:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "G4HEPEvtParticle.hh"

G4Allocator<G4HEPEvtParticle> aHEPEvtParticleAllocator;

G4HEPEvtParticle::G4HEPEvtParticle()
{;}

G4HEPEvtParticle::G4HEPEvtParticle(G4PrimaryParticle* pp,
        G4int isthep, G4int jdahep1, G4int jdahep2)
:theParticle(pp),ISTHEP(isthep),JDAHEP1(jdahep1),JDAHEP2(jdahep2)
{;}

G4HEPEvtParticle::~G4HEPEvtParticle()
{;}

const G4HEPEvtParticle & 
G4HEPEvtParticle::operator=(const G4HEPEvtParticle &right)
{ return *this; }

int G4HEPEvtParticle::operator==(const G4HEPEvtParticle &right) const
{ return false; }
int G4HEPEvtParticle::operator!=(const G4HEPEvtParticle &right) const
{ return true; }




