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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2SimpleReporter.cc,v 1.3 2001-07-11 10:02:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#include "tst2SimpleReporter.hh"
#include "G4ios.hh"
#include "globals.hh"

 tst2SimpleReporter::tst2SimpleReporter():tst2VParticleReporter()
{
 
}

 tst2SimpleReporter::~tst2SimpleReporter()
{
}    

 void tst2SimpleReporter::Print(const tst2ParticleContainer& container, 
                               const G4String& option)
{
  pList = &container;
  G4cout << " Encoding    " << "name " << G4endl;
  for (G4int i=0; i< entries(); i++){
	G4cout << GetEncoding(i) << " :   " << GetParticle(i)->GetParticleName() << G4endl;
  }
}    
