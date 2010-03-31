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
// $Id: T07PhysicsList.hh,v 1.7 2010-03-31 06:51:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07PhysicsList_h
#define T07PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4PhotoElectricEffect;
class G4ComptonScattering;
class G4GammaConversion;

class G4eMultipleScattering;
class G4MuMultipleScattering;
class G4hMultipleScattering;

class G4eIonisation;
class G4eBremsstrahlung;
class G4eplusAnnihilation;

class G4MuIonisation;
class G4MuBremsstrahlung;
class G4MuPairProduction;

class G4hIonisation;

class T07PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07PhysicsList: public G4VUserPhysicsList
{
  public:
    T07PhysicsList();
    virtual ~T07PhysicsList();

  protected:
    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
 
    virtual void SetCuts();
    
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
    
  public:
  // this method allow to set on/off the processes
    void SetStatusEmProcess();
    
  private:

    G4PhotoElectricEffect* thePhotoElectricEffect;
    G4ComptonScattering*   theComptonScattering;
    G4GammaConversion*     theGammaConversion;
    
    G4eMultipleScattering* theeminusMultipleScattering;
    G4eIonisation*         theeminusIonisation;
    G4eBremsstrahlung*     theeminusBremsstrahlung;
    
    G4eMultipleScattering* theeplusMultipleScattering;
    G4eIonisation*         theeplusIonisation;
    G4eBremsstrahlung*     theeplusBremsstrahlung;
    G4eplusAnnihilation*   theeplusAnnihilation;
    
    T07PhysicsListMessenger* physicsListMessenger;
};

#endif
















