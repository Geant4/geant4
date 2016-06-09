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
//  Code developed by: Susanna Guatelli
//
// $Id: BrachyFactory.hh,v 1.1 2004/05/25 07:32:35 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//    **********************************
//    *                                *
//    *      BrachyFactory.hh          *
//    *                                *
//    **********************************
//
// this is the abstract class which manages the sources' realisation
// in terms of geometry and primary particle generated in the radioactive core

#ifndef BrachyFactory_h
#define BrachyFactory_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class G4VPhysicalVolume;
class BrachyAnalysisManager;

class BrachyFactory
{
public:

  BrachyFactory();
  virtual ~BrachyFactory();

  //Source primary particles' management ...
  virtual G4VUserPrimaryGeneratorAction* CreatePrimaryGeneratorAction() = 0;

  // ...
  // this function manages the creation of the source in terms of 
  // geometry ...
  virtual void CreateSource(G4VPhysicalVolume*) = 0;

  // ...
  // This function deletes the source in the geometry (phantom); to be
  // used when the user needs to substitute one source with another one ...
  virtual void CleanSource() = 0;
};
#endif
