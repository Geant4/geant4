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
// $Id: G4RToEConvForProton.cc,v 1.2 2002-12-16 11:15:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    5 Oct. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4RToEConvForProton.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "g4std/iomanip"
#include "g4std/strstream"

G4RToEConvForProton::G4RToEConvForProton() : G4VRangeToEnergyConverter()
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("proton");
  if (theParticle ==0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << " G4RToEConvForProton::G4RToEConvForProton() ";
      G4cout << " proton is not defined !!" << G4endl;
    }
#endif
  } 
}

G4RToEConvForProton::~G4RToEConvForProton()
{ 
}

