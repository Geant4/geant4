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
// ClassName:   G4HyperonFTFPBuilder
//
// Author: 2012 G.Folger
//
// Modified:
// 06.05.2020 A.Ribon : introduced G4VHyperonBuilder
// 12.04.2017 A.Dotti : move to new design with base class
//
//----------------------------------------------------------------------------

#ifndef G4HyperonFTFPBuilder_h
#define G4HyperonFTFPBuilder_h 1

#include "G4VHyperonBuilder.hh"
#include "globals.hh"

class G4TheoFSGenerator;
class G4CascadeInterface;
class G4VCrossSectionDataSet;


class G4HyperonFTFPBuilder : public G4VHyperonBuilder {
  public: 
    G4HyperonFTFPBuilder( G4bool quasiElastic = false );
    virtual ~G4HyperonFTFPBuilder();

    virtual void Build( G4HadronElasticProcess* ) final override {}
    virtual void Build( G4HadronInelasticProcess* aP ) final override;

    // The energy limits refer to the string model FTF:
    // -  the max energy is the same for hyperons and hyperons;
    // -  the min energy is for hyperons only (0.0 is assumed for antihyperons)
    virtual void SetMinEnergy( G4double val ) final override { theMin = val; }
    virtual void SetMaxEnergy( G4double val ) final override { theMax = val; }

    using G4VHyperonBuilder::Build;  // Prevent compiler warning

  private: 
    G4TheoFSGenerator* theHyperonFTFP;
    G4TheoFSGenerator* theAntiHyperonFTFP;
    G4CascadeInterface* theBertini;
    G4VCrossSectionDataSet* theInelasticCrossSection;
    G4double theMin;  // Min energy for FTF for hyperons only (0.0 is assumed for antihyperons)
    G4double theMax;  // Max energy for FTF for hyperons and antihyperons
};

#endif
