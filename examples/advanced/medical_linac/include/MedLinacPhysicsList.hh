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
//// $Id: MedLinacPhysicsList.hh,v 1.4 2004/06/18 09:17:39 gunter Exp $
//
//
// Code developed by: M. Piergentili



#ifndef MedLinacPhysicsList_h
#define MedLinacPhysicsList_h 1

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

  //private:
    
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







