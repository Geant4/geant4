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
// $Id: G4HadronPhysicsFTF_BIC.hh 105736 2017-08-16 13:01:11Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTF_BIC_h
#define G4HadronPhysicsFTF_BIC_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"


#include "G4Cache.hh"

class G4ComponentGGHadronNucleusXsc;
class G4VCrossSectionDataSet;

class G4HadronPhysicsFTF_BIC : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsFTF_BIC(G4int verbose =1);
    G4HadronPhysicsFTF_BIC(const G4String& name,G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsFTF_BIC();

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    virtual void TerminateWorker() override;
  protected:
    G4bool QuasiElastic;
    //This calls the specific ones for the different particles in order
    virtual void CreateModels();
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();
    virtual void Kaon();
    virtual void Others();
    virtual void DumpBanner() {}
    //This contains extra configurataion specific to this PL
    virtual void ExtraConfiguration();
    
    G4double maxBIC_pion;
    G4double maxBERT_kaon; //Bertini for kaons
    G4double maxBIC_proton;
    G4double maxBIC_neutron;
   
    //Thread-private data
    G4VectorCache<G4VCrossSectionDataSet*> xs_ds;
    G4Cache<G4ComponentGGHadronNucleusXsc*> xs_k;
};

#endif

