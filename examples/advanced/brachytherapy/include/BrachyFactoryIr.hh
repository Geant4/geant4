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
// $Id: BrachyFactoryIr.hh,v 1.2 2002-11-18 15:18:36 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include"BrachyFactory.hh"
#include "G4RunManager.hh"
class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyFactory;
class BrachyPrimaryGeneratorActionIr;
class BrachyDetectorConstructionIr;

class BrachyFactoryIr:public BrachyFactory
{
public:
  BrachyFactoryIr();
 ~BrachyFactoryIr();
 G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction();
  void CreateSource(G4VPhysicalVolume*);
  void CleanSource();
private:
  BrachyDetectorConstructionIr* pIridio;
};
#endif
