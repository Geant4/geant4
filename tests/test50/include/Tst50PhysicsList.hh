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
// $Id: Tst50PhysicsList.hh,v 1.6 2003-02-07 13:27:49 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst50PhysicsList_h
#define Tst50PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;
class Tst50PhysicsListMessenger;


class Tst50PhysicsList: public G4VUserPhysicsList
{
  public:
    Tst50PhysicsList(G4bool,G4bool,G4bool,G4bool,G4bool);
   ~Tst50PhysicsList();
public:
  void SetRangeConditions(G4String );
  void SetElectronCut(G4double);

 protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();

   
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
  void ConstructHad();
private: 
  G4LowEnergyIonisation*  loweIon;
  G4LowEnergyPhotoElectric* lowePhot;
  G4LowEnergyBremsstrahlung* loweBrem;
  Tst50PhysicsListMessenger* physicsListMessenger;
  G4String rangeOn; 
   G4double cutForElectron;
  G4bool Low;
  G4bool  RangeOn;
  G4bool Stopping;
  G4bool RadiationY;
  G4bool Adronic;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

 





