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
// $Id: BrachyFactoryIr.hh 69765 2013-05-14 10:11:22Z gcosmo $
//
//    **********************************
//    *                                *
//    *      BrachyFactoryIr.hh        *
//    *                                *
//    **********************************
//
//
//
#ifndef BrachyFactoryIr_h
#define BrachyFactoryIr_h 1

#include "BrachyFactory.hh"
#include "G4RunManager.hh"

class BrachyFactory;
class BrachyDetectorConstructionFlexi;
class G4ParticleGun;
class G4Run;
class G4Event;

// This class manages the creation of the Nucletron Flexi Source

// This class manages the creation of iridum source used in endocavitary
// brachytherapy ...
class BrachyFactoryIr : public BrachyFactory
{
public:
  explicit BrachyFactoryFlexi();
  ~BrachyFactoryFlexi();
  void CreateSource(G4VPhysicalVolume*)override;
  void CleanSource() override;

private:
  BrachyDetectorConstructionFlexi* fFlexiSource; 
};
#endif
