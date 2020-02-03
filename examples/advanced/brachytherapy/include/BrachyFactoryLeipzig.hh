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
//
//    **********************************
//    *                                *
//    *      BrachyFactoryLeipzig.hh   *
//    *                                *
//    **********************************
//code developed by: Susanna Guatelli
//
// This class manages the creation of iridum source used in superficial
// brachytherapy ...

#ifndef BrachyFactoryLeipzig_h
#define BrachyFactoryLeipzig_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "BrachyFactory.hh"
#include "G4RunManager.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyFactory;
class BrachyDetectorConstructionLeipzig;

class BrachyFactoryLeipzig : public BrachyFactory
{
public:
  BrachyFactoryLeipzig();
 ~BrachyFactoryLeipzig();

  void CreateSource(G4VPhysicalVolume*);
  void CleanSource();

private:
  BrachyDetectorConstructionLeipzig* leipzigSource;
};
#endif
