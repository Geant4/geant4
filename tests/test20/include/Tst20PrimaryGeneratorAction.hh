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
// $Id: Tst20PrimaryGeneratorAction.hh,v 1.5 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#ifndef Tst20PrimaryGeneratorAction_h
#define Tst20PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Tst20DetectorConstruction;
//class Tst20PrimaryGeneratorMessenger;


class Tst20PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  Tst20PrimaryGeneratorAction(Tst20DetectorConstruction*);    
  ~Tst20PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String value) { rndmFlag = value;}
  void SetXvertex(G4double x) ;
  void SetYvertex(G4double y) ;
  void SetZvertex(G4double z) ;

  G4String GetPrimaryName();                

private:

  G4ParticleGun* particleGun;	//pointer a to G4 service class
  Tst20DetectorConstruction* detector; //pointer to the geometry
    
  //  Tst20PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  G4String rndmFlag;	//flag for a random impact point       

  G4String primaryParticleName ;
  G4double xVertex;
  G4double yVertex;
  G4double zVertex;
  G4bool vertexDefined ;

};

#endif


