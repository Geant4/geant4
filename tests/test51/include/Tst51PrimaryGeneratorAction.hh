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
// $Id: Tst51PrimaryGeneratorAction.hh,v 1.1 2005-07-05 11:05:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
 
#ifndef Tst51PrimaryGeneratorAction_h
#define Tst51PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4ParticleGun;
class G4Event;
//class Tst51PrimaryGeneratorMessenger; 

class Tst51PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst51PrimaryGeneratorAction();    
  ~Tst51PrimaryGeneratorAction();
  void GeneratePrimaries(G4Event*);
  
private:
  G4ParticleGun* particleGun;
  //Tst51PrimaryGeneratorMessenger* gunMessenger; 
};
#endif










