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
// $Id: Tst50Particles.cc,v 1.3 2004-06-02 09:46:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Created
//
// -------------------------------------------------------------------

#include "Tst50Particles.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4IonO.hh"
#include "G4IonC12.hh"
#include "G4IonSi28.hh"
#include "G4IonFe52.hh"
Tst50Particles::Tst50Particles(const G4String& name)
  :  G4VPhysicsConstructor(name)
{ }

Tst50Particles::~Tst50Particles()
{}

void Tst50Particles::ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Proton :: ProtonDefinition(); 
  G4Alpha:: AlphaDefinition(); 
  G4IonO::IonODefinition();  
  G4IonC12::IonC12Definition();
  G4IonSi28::IonSi28Definition();
  G4IonFe52::IonFe52Definition();
}
