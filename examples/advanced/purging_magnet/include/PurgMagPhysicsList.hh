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
// Code developed by:
//  S.Larsson
//
//    ********************************
//    *                              *
//    *    PurgMagPhysicsList.hh     *
//    *                              *
//    ********************************
//
//
// $Id: PurgMagPhysicsList.hh,v 1.2 2004/06/18 09:17:46 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagPhysicsList_h
#define PurgMagPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4GammaConversion;     
class G4ComptonScattering;   
class G4PhotoElectricEffect;       

class G4eIonisation;          
class G4eBremsstrahlung;     


class PurgMagPhysicsList: public G4VUserPhysicsList
{
public:
  PurgMagPhysicsList();
  ~PurgMagPhysicsList();
  
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
  void SetProtonCut(G4double);
  
  void SetGammaLowLimit(G4double);
  void SetElectronLowLimit(G4double);
  void SetPositronLowLimit(G4double);
  void SetProtonLowLimit(G4double);
  void SetGEPLowLimit(G4double);

  void SetGELowLimit(G4double);
  
private:
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructBarions();
  
protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  
private:
  
};
#endif



