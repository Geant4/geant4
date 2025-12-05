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
// ClassName:   G4FTFPProtonBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 30.03.2009 V.Ivanchenko create cross section by new
//
//----------------------------------------------------------------------------
//
#ifndef G4FTFPProtonBuilder_h
#define G4FTFPProtonBuilder_h 

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4VProtonBuilder.hh"

#include "G4TheoFSGenerator.hh"

class G4FTFPProtonBuilder : public G4VProtonBuilder
{
  public: 
    G4FTFPProtonBuilder(G4bool quasiElastic=false);
    ~G4FTFPProtonBuilder() override = default;

    void Build(G4HadronElasticProcess *) override {};
    void Build(G4HadronInelasticProcess * aP) override;
    
    void SetMinEnergy(G4double aM) override { theMin = aM; };
    void SetMaxEnergy(G4double aM) override { theMax = aM; };

    using G4VProtonBuilder::Build; //Prevent compiler warning

  private:
    G4TheoFSGenerator * theModel;
    G4double theMin;
    G4double theMax;

};

#endif

