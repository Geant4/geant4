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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////
#ifndef IORTPrimaryGeneratorAction_h
#define IORTPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class IORTPrimaryGeneratorMessenger;
class IORTPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  IORTPrimaryGeneratorAction();    
  ~IORTPrimaryGeneratorAction();
  
public:
  // Methods to change the parameters of primary particle generation 
  // interactively
  void SetsigmaEnergy(G4double);
  void SetmeanKineticEnergy(G4double);
  void GeneratePrimaries(G4Event*);
  void SetXposition(G4double);
  void SetYposition(G4double);
  void SetZposition(G4double);
  void SetsigmaY(G4double);
  void SetsigmaZ(G4double);
  // void SetsigmaMomentumY(G4double);
  // void SetsigmaMomentumZ(G4double);
  void SetTheta(G4double); // aggiunto
  G4double GetmeanKineticEnergy(void);
  G4ParticleGun *GetParticleGun(void){return particleGun;}
    
private:
  void SetDefaultPrimaryParticle();
  G4double meanKineticEnergy;
  G4double sigmaEnergy;
  G4double X0;
  G4double Y0;
  G4double Z0;
  G4double sigmaY;
  G4double sigmaZ;
  //  G4double sigmaMomentumY;
  //  G4double sigmaMomentumZ;
  G4double Theta;  // aggiunto

private:
  G4ParticleGun*    		          particleGun;
  IORTPrimaryGeneratorMessenger* gunMessenger; 
};

#endif


