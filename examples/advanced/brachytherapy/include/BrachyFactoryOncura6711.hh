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
//    *   BrachyFactoryOncura6711.hh   *
//    *                                *
//    **********************************
//
// Author: D. Cutajar, A. Le
//
#ifndef BrachyFactoryOncura6711_h
#define BrachyFactoryOncura6711_h 1

#include "BrachyFactory.hh"
#include "G4RunManager.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyFactory;
class BrachyDetectorConstructionOncura6711;

// This class manages the creation of iodine source
// brachytherapy ...
class BrachyFactoryOncura6711 : public BrachyFactory
{
public:
  BrachyFactoryOncura6711();
  ~BrachyFactoryOncura6711();
  void CreateSource(G4VPhysicalVolume*);
  void CleanSource();

private:
  BrachyDetectorConstructionOncura6711* Oncura6711IodineSource;
};
#endif
