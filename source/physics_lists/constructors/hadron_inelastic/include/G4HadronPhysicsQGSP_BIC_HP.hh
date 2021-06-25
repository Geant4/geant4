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
// ClassName:   G4HadronPhysicsQGSP_BIC_HP
//
// Author: 2006 G.Folger
//
// Based on G4HadronPhysicsQGSP_BIC
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 25.04.2007 G.Folger: Add quasielastic as option, use quasielastic by default
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
// 05.05.2020 A.Ribon: Use QGSP for antibaryons at high energies
// 07.05.2020 A.Ribon: Use QGSP for hyperons (and anti-hyperons) at high energies
// 20.05.2020 A.Ribon: Refactoring of the class (keeping same functionalities)
// 02.10.2020 V.Ivanchenko: added c-,b- particles and cross section biasing
//
//----------------------------------------------------------------------------

#ifndef G4HadronPhysicsQGSP_BIC_HP_h
#define G4HadronPhysicsQGSP_BIC_HP_h 1

#include "G4HadronPhysicsQGSP_BIC.hh"


class G4HadronPhysicsQGSP_BIC_HP : public G4HadronPhysicsQGSP_BIC {
  public: 
    G4HadronPhysicsQGSP_BIC_HP( G4int verbose = 1 );
    G4HadronPhysicsQGSP_BIC_HP( const G4String& name, G4bool quasiElastic = true );
    virtual ~G4HadronPhysicsQGSP_BIC_HP() {};

    // copy constructor and hide assignment operator
    G4HadronPhysicsQGSP_BIC_HP(G4HadronPhysicsQGSP_BIC_HP &) = delete;
    G4HadronPhysicsQGSP_BIC_HP & operator =
    (const G4HadronPhysicsQGSP_BIC_HP &right) = delete;

  protected:
    void Neutron() override;

};

#endif
