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
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from G4HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
// 18.10.2020 V.Ivanchenko use inheritance from G4HadronPhysicsQGSP_BERT
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsQGS_BIC_h
#define G4HadronPhysicsQGS_BIC_h 1

#include "G4HadronPhysicsQGSP_BERT.hh"

class G4HadronPhysicsQGS_BIC : public G4HadronPhysicsQGSP_BERT
{
  public: 
    G4HadronPhysicsQGS_BIC(G4int verbose =1);
    G4HadronPhysicsQGS_BIC(const G4String& name, G4bool quasiElastic=true);
    virtual ~G4HadronPhysicsQGS_BIC();

    // copy constructor and hide assignment operator
    G4HadronPhysicsQGS_BIC(G4HadronPhysicsQGS_BIC &) = delete;
    G4HadronPhysicsQGS_BIC & operator =
    (const G4HadronPhysicsQGS_BIC &right) = delete;

  protected:
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();

  private:
    G4double minBERT_pion;
    G4double maxBIC_pion;
    
};

#endif

