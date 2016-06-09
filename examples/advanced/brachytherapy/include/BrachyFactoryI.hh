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
//    *******************************
//    *                             *
//    *    BrachyFactoryI.cc        *
//    *                             *
//    *******************************
//
//
// Code developed by:
//  S.Guatelli
//
// 
// $Id: BrachyFactoryI.hh,v 1.5 2003/05/22 17:20:41 guatelli Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// 
// --------------------------------------------------------------
#ifndef BrachyFactoryI_h
#define BrachyFactoryI_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "BrachyDetectorConstructionI.hh"
#include "BrachyFactory.hh"
#include "G4RunManager.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyFactory;
class BrachyPrimaryGeneratorActionI;
class BrachyDetectorConstructionI;

// This class manages the creation of Bebig Isoseed I-125 source 
// used in interstitial brachytherapy ...
class BrachyFactoryI:public BrachyFactory
{
public:
  BrachyFactoryI();
 ~BrachyFactoryI();

  G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction();
  void CreateSource(G4VPhysicalVolume*);
  void CleanSource();

private:
  BrachyDetectorConstructionI* iodiumSource;
};
#endif
