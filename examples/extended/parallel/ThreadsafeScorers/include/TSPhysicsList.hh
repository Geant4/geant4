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
/// \file parallel/ThreadsafeScorers/include/TSPhysicsList.hh
/// \brief Definition of the TSPhysicsList class
//
//
// $Id: TSPhysicsList.hh 93110 2015-11-05 08:37:42Z jmadsen $
//
//
/// This is a very, very extensive physics list and step-limiters are applied
///     to many particles. The reasoning behind this is because we wan't to put
///     as much pressure on the atomics as possible and produce as much
///     round-off error as possible. See descriptions in README and
///     TSDetectorConstruction for more details.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef tsphysicslist_hh
#define tsphysicslist_hh 1


#include "globals.hh"
#include "G4VUserPhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"

class TSPhysicsList : public G4VUserPhysicsList
{
public:
    typedef std::deque<G4VPhysicsConstructor*> PhysicsSet_t;

public:
    TSPhysicsList();
    virtual ~TSPhysicsList();

    static TSPhysicsList* Instance();

public:
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

private:
    static TSPhysicsList* fgInstance;

    G4EmStandardPhysics_option4*    fEmPhysics_opt4;
    G4DecayPhysics*                 fDecayPhysics;
    G4RadioactiveDecayPhysics*      fRadDecayPhysics;
    G4HadronPhysicsQGSP_BERT_HP*    fHadronInelasticPhysics;
    G4HadronElasticPhysicsHP*       fHadronElasticPhysics;
    G4IonElasticPhysics*            fIonElasticPhysics;
    G4IonBinaryCascadePhysics*      fIonBinaryCascadePhysics;

    PhysicsSet_t fConstructors;

    G4double fDefaultCutValue;

};

#endif
