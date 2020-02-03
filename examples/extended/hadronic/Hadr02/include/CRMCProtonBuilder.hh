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
/// \file hadronic/Hadr02/include/CRMCProtonBuilder.hh
/// \brief Definition of the CRMCProtonBuilder class
//
//
//---------------------------------------------------------------------------
//
// ClassName: CRMCProtonBuilder
//
// Author:    2018 Alberto Ribon
//
// Physics builder to treat the final state of inelastic proton-nuclear
// interactions with the wrapper hadronic model around CRMC.
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef CRMCProtonBuilder_h
#define CRMCProtonBuilder_h 

#include "globals.hh"
#include "G4VProtonBuilder.hh"
#include "G4CRMCModel.hh"


class CRMCProtonBuilder : public G4VProtonBuilder {
  public:
    CRMCProtonBuilder();
    virtual ~CRMCProtonBuilder();
    virtual void Build( G4HadronElasticProcess* aP ) final override;
    virtual void Build( G4ProtonInelasticProcess* aP ) final override;
    inline void SetMinEnergy( G4double aM ) final override { fMin = aM; }
    inline void SetMaxEnergy( G4double aM ) final override { fMax = aM; }
    using G4VProtonBuilder::Build;  // Prevent compiler warning
  private:
    G4CRMCModel* fModel;    
    G4double fMin;
    G4double fMax;
};

#endif

