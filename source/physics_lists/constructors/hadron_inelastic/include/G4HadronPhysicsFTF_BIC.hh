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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
// 19.06.2008 G.Folger: change default for QE to NOT use Chips QE
// 18.10.2020 V.Ivanchenko use inheritance from G4HadronPhysicsFTFP_BERT
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTF_BIC_h
#define G4HadronPhysicsFTF_BIC_h 1

#include "G4HadronPhysicsFTFP_BERT.hh"

class G4HadronPhysicsFTF_BIC : public G4HadronPhysicsFTFP_BERT
{
  public: 
    G4HadronPhysicsFTF_BIC(G4int verbose =1);
    G4HadronPhysicsFTF_BIC(const G4String& name, G4bool quasiElastic=false);
    ~G4HadronPhysicsFTF_BIC() override;

    // copy constructor and hide assignment operator
    G4HadronPhysicsFTF_BIC(G4HadronPhysicsFTF_BIC &) = delete;
    G4HadronPhysicsFTF_BIC & operator =
    (const G4HadronPhysicsFTF_BIC &right) = delete;

  protected:
    void Neutron() override;
    void Proton() override;
    void Pion() override;
    void Kaon() override;

  private:
    G4double maxBIC_pion;
    G4double minBERT_pion;
};

#endif

