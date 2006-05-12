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
// $Id: BrachyPrimaryGeneratorAction.hh,v 1.11 2006-05-12 17:08:06 guatelli Exp $
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


