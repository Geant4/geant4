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
// $Id: NTSTPhysicsList.hh,v 1.2 2003-12-09 15:35:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTPhysicsList_h
#define NTSTPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class NTSTLooperDeath;

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4MultipleScattering;

class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4MuIonisation;
class G4MuBremsstrahlung;
class G4MuPairProduction;

class G4hIonisation;

class NTSTPhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTPhysicsList: public G4VUserPhysicsList
{
public:
  NTSTPhysicsList();
  ~NTSTPhysicsList();
  //
  // set methods
  //
  inline void SetBgsTran(const G4bool on) { useBgsTran = on; }
  void SetMinimumEnergyCut(const G4double e);
  void SetMaximumEnergyCut(const G4double e);
  inline void SetLengthCut(const G4double e) { Cut = e ; }
  void SetLooperCut(const G4double e);
    
protected:
  void AddBgsTransportation();
  
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
 
  void SetCuts();
    
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBarions();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructLeptHad();
  void ConstructHad();
    
public:
  // this method allow to set on/off the processes
  void SetStatusEmProcess();
    
private:
  G4bool   useBgsTran;
  G4double MinimumEnergyCut;
  G4double MaximumEnergyCut;
  G4double Cut;
  G4double LooperCut;

  // processes

  NTSTLooperDeath*       theLooperDeath;

  G4PhotoElectricEffect* thePhotoElectricEffect;
  G4ComptonScattering*   theComptonScattering;
  G4GammaConversion*     theGammaConversion;
    
  G4MultipleScattering*  theeminusMultipleScattering;
  G4eIonisation*         theeminusIonisation;
  G4eBremsstrahlung*     theeminusBremsstrahlung;
    
  G4MultipleScattering*  theeplusMultipleScattering;
  G4eIonisation*         theeplusIonisation;
  G4eBremsstrahlung*     theeplusBremsstrahlung;
  G4eplusAnnihilation*   theeplusAnnihilation;
    
  NTSTPhysicsListMessenger* physicsListMessenger;
};

#endif



