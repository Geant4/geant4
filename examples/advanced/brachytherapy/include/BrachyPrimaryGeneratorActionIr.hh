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
// $Id: BrachyPrimaryGeneratorActionIr.hh,v 1.2 2002-11-18 15:18:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    *************************************************
//    *                                               *
//    *     BrachyPrimaryGeneratorActionIr.hh"        *
//    *                                               *
//    *************************************************
//
// 
//Primary particles of the Iridium source. They are delivered from a random
//point of the radionuclide  with random direction. 

#ifndef BrachyPrimaryGeneratorActionIr_h
#define BrachyPrimaryGeneratorActionIr_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RadioactiveDecay.hh"
#include "BrachyPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyPrimaryGeneratorAction;
class BrachyPrimaryGeneratorActionIr: public G4VUserPrimaryGeneratorAction
{
 public:
	BrachyPrimaryGeneratorActionIr();
      ~BrachyPrimaryGeneratorActionIr();

 public:
      void GeneratePrimaries(G4Event* anEvent);
  G4double GetEnergy(){return Energy;};
     
 private:
  // G4ParticleGun* m_pParticleGun;
     G4ParticleGun* m_pParticleGun;
  //G4RadioactiveDecay *m_pRadioactiveDecay;
        G4double Energy;
       G4RadioactiveDecay *m_pRadioactiveDecay;
  //    BrachyPrimaryGeneratorMessenger* gunMessenger; 
  G4std::vector<G4double> vettore;
};

#endif






























































