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
// $Id: RemSimVPrimaryGeneratorFactory.hh,v 1.8 2004/11/22 16:51:38 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#ifndef RemSimVPrimaryGeneratorFactory_h
#define RemSimVPrimaryGeneratorFactory_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class RemSimVPrimaryGeneratorFactory : public G4VUserPrimaryGeneratorAction
{
public:
  RemSimVPrimaryGeneratorFactory();
  virtual ~RemSimVPrimaryGeneratorFactory();

public:
  virtual void GeneratePrimaries(G4Event* anEvent)=0;
  virtual G4double GetInitialEnergy()=0;
  virtual void SetMoon(G4bool)=0; // Set moon configuration 
  virtual void SetParticle(G4String)=0;
};
#endif


