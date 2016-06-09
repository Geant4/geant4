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
//// $Id: MedLinacPhysicsList.hh,v 1.7 2006/06/29 16:03:57 gunter Exp $
//
//
// Code developed by: M. Piergentili



#ifndef MedLinacPhysicsList_h
#define MedLinacPhysicsList_h 1
#include "MedLinacPhysicsListMessenger.hh"

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4ProductionCuts;
class G4LowEnergyPhotoElectric;
class G4LowEnergyIonisation;
class G4LowEnergyBremsstrahlung;

class MedLinacPhysicsList: public G4VUserPhysicsList
{
  public:
    MedLinacPhysicsList();
    ~MedLinacPhysicsList();

  protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

public:
  void SetCut (G4double);
  G4double GetCut()  {return defaultCut;};
  
private:
  G4double defaultCut;
  MedLinacPhysicsListMessenger* physicsListMessenger;

  //G4double cutForGamma;
  //G4double cutForElectron;
  //G4double cutForPositron;

protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
 
  protected:
  // these methods Construct physics processes and register them
  //void ConstructGeneral();
  void ConstructEM();

  // private:

  G4LowEnergyIonisation*  loweIon;
  G4LowEnergyPhotoElectric* lowePhot;
  G4LowEnergyBremsstrahlung* loweBrem;
};

#endif







