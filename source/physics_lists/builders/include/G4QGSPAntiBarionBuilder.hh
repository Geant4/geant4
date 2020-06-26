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
//--------------------------------------------------------------------------
// ClassName: G4QGSPAntiBarionBuilder
// Author: Alberto Ribon
// Date: May 2020
// Description: builder class that uses QGSP for anti_proton and
//              anti_neutron, while, for the time being, delegates FTFP for
//              the light anti-ions (anti_deuteron, anti_triton, anti_He3
//              and anti_alpha).
//              For safety, better to use this class above >~ 12 GeV.
// Modified:
//---------------------------------------------------------------------------

#ifndef G4QGSPAntiBarionBuilder_h
#define G4QGSPAntiBarionBuilder_h 1

#include "globals.hh"
#include "G4VAntiBarionBuilder.hh"

class G4HadronElasticProcess;
class G4VCrossSectionDataSet;
class G4TheoFSGenerator;


class G4QGSPAntiBarionBuilder : public G4VAntiBarionBuilder {
  public: 
    G4QGSPAntiBarionBuilder( G4bool quasiElastic = false );
    virtual ~G4QGSPAntiBarionBuilder() {};

    virtual void Build( G4HadronElasticProcess* ) final override {}
    virtual void Build( G4AntiProtonInelasticProcess* aP ) final override;
    virtual void Build( G4AntiNeutronInelasticProcess* aP ) final override;
    virtual void Build( G4AntiDeuteronInelasticProcess* ) final override;
    virtual void Build( G4AntiTritonInelasticProcess* ) final override;
    virtual void Build( G4AntiHe3InelasticProcess* ) final override;
    virtual void Build( G4AntiAlphaInelasticProcess* ) final override;
    
    virtual void SetMinEnergy( G4double val ) final override { theMin = val; }
    virtual void SetMaxEnergy( G4double val ) final override { theMax = val; }

    using G4VAntiBarionBuilder::Build;  // Prevent compiler warning

  private:
    G4TheoFSGenerator* theQGSmodel;
    G4TheoFSGenerator* theFTFmodel;
    G4VCrossSectionDataSet* theAntiNucleonData;
    G4double theMin;
    G4double theMax;
};

#endif
