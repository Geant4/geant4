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
// $Id: BrachyPrimaryGeneratorAction.hh,v 1.12 2006-06-29 15:47:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    ********************************************
//    *                                          *
//    *      BrachyPrimaryGeneratorAction.hh     *
//    *                                          *
//    ********************************************

// This class must be implemented because it is mandatory

#ifndef BrachyPrimaryGeneratorAction_h
#define BrachyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyFactory;
class  BrachyPrimaryGeneratorMessenger;

class BrachyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
   BrachyPrimaryGeneratorAction();
   ~BrachyPrimaryGeneratorAction();

 public:
  void GeneratePrimaries(G4Event* anEvent);
  void SwitchEnergy(G4String);

private:
  BrachyFactory* factory;
  BrachyPrimaryGeneratorMessenger* primaryMessenger;
};

#endif


