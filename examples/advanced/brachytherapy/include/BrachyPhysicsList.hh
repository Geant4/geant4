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
// $Id: BrachyPhysicsList.hh,v 1.6 2002-11-18 15:18:36 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    **********************************
//    *                                *
//    *      BrachyPhysicsList.hh      *
//    *                                *
//    **********************************

#ifndef BrachyPhysicsList_h
#define BrachyPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

class BrachyPhysicsList: public G4VUserPhysicsList
{
  public:
    BrachyPhysicsList();
   ~BrachyPhysicsList();

  protected:
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
 
    void SetCuts();
  
  public: 
    // Set Cuts
    void SetGammaCut(G4double);
    void SetElectronCut(G4double);
    void SetPositronCut(G4double);
    
    void SetGammaLowLimit(G4double);
    void SetElectronLowLimit(G4double);
    void SetGELowLimit(G4double);
    void SetLowEnSecPhotCut(G4double);
    void SetLowEnSecElecCut(G4double);
    
  private:
    
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;
    
  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    
  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();

  private:
  G4LowEnergyIonisation*  loweIon;
  G4LowEnergyPhotoElectric* lowePhot;
  G4LowEnergyBremsstrahlung* loweBrem;
  
};

#endif



