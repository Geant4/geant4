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
// $Id: G4GammaNuclearReaction.hh,v 1.15 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4GammaNuclearReaction - header file for CHIPS
// Created: J.P. Wellisch, following M.Kossov's algorithm. 2000/08/18 
// The last update: J.P. Wellisch, Thu Jun  6 2002.
// 01.09.2008 V.Ivanchenko move inline to source and define interaction name
// 17.02.2009 M.Kossov, now it is recommended to use the G4QCollision process

#ifndef G4GammaNuclearReaction_h
#define G4GammaNuclearReaction_h 1

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4GammaNuclearReaction : public G4HadronicInteraction
{
public: 

  G4GammaNuclearReaction();

  virtual ~G4GammaNuclearReaction();
    
  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                         G4Nucleus& aTargetNucleus);

private:

  G4ChiralInvariantPhaseSpace theModel;

};

#endif
