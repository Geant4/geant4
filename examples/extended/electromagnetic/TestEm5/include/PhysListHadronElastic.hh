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
// $Id: PhysListHadronElastic.hh,v 1.3 2004/06/21 10:57:11 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Class Description:
//      This class is an derived class of G4VPhysicsConstructor
//      It is provide PhysicsList for hadron eleastic process
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysListHadronElastic_h
#define PhysListHadronElastic_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4HadronElasticProcess.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysListHadronElastic : public G4VPhysicsConstructor
{
  public: 
    PhysListHadronElastic(const G4String& name = "elastic");
    virtual ~PhysListHadronElastic();

  public: 
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    void ConstructParticle() {};
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    void ConstructProcess();

  private:
    // Elastic Process
    G4HadronElasticProcess theElasticProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif








