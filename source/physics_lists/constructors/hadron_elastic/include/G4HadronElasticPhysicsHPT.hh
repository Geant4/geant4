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
// ClassName:   G4HadronElasticPhysicsHPT
//
// Author: Alberto Ribon (CERN), April 2023
//
// This class inherits from G4HadronElasticPhysicsHP and activates
// the special treatment of elastic scattering of thermal neutrons
// (i.e. with kinetic energy below 4 eV).
// This special treatment, called Thermal Scattering Law (TSL), is based on
// the S(alpha, beta) approach, which relies on both experimental measurements
// and molecular dynamics calculations.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronElasticPhysicsHPT_h
#define G4HadronElasticPhysicsHPT_h 1

#include "G4HadronElasticPhysicsHP.hh"


class G4HadronElasticPhysicsHPT : public G4HadronElasticPhysicsHP {
  public: 
    explicit G4HadronElasticPhysicsHPT( G4int ver = 1 ); 
    virtual ~G4HadronElasticPhysicsHPT() = default;
    void ConstructProcess() final;
    G4HadronElasticPhysicsHPT( G4HadronElasticPhysicsHPT & ) = delete;
    G4HadronElasticPhysicsHPT & operator=( const G4HadronElasticPhysicsHPT &right ) = delete;
};


#endif
