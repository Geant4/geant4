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
/// \file hadronic/Hadr02/include/HadronPhysicsCRMC_FTFP_BERT.hh
/// \brief Definition of the HadronPhysicsCRMC_FTFP_BERT class
//
//
//---------------------------------------------------------------------------
//
// ClassName: HadronPhysicsCRMC_FTFP_BERT
//
// Author:    2018 Alberto Ribon
//
// This is a variant of HadronPhysicsFTFP_BERT whereby CRMC is used 
// for modeling final-state for pion- , kaon- , proton- and neutron-
// nuclear inelastic interactions at very high energies.
// For other hadron projectile types (e.g. hyperons, antinucleons and
// antihyperons) the usual FTFP_BERT approach is used at all energies.
// The inelastic hadronic cross sections are, for all hadron projectiles
// and energies, the usual ones (exactly as in FTFP_BERT).
//
// Modified:
// -  18-May-2021 Alberto Ribon : Migrated to newer physics constructor
//                                and used the latest Geant4-CRMC interface.
//
//----------------------------------------------------------------------------
//
#ifndef HadronPhysicsCRMC_FTFP_BERT_h
#define HadronPhysicsCRMC_FTFP_BERT_h 1

#include "G4HadronPhysicsFTFP_BERT.hh"


class HadronPhysicsCRMC_FTFP_BERT : public G4HadronPhysicsFTFP_BERT {
  public: 
    HadronPhysicsCRMC_FTFP_BERT( G4int verbose = 1 );
    HadronPhysicsCRMC_FTFP_BERT( const G4String& name, G4bool quasiElastic = false );
    ~HadronPhysicsCRMC_FTFP_BERT() override;

    // copy constructor and hide assignment operator
    HadronPhysicsCRMC_FTFP_BERT( HadronPhysicsCRMC_FTFP_BERT & ) = delete;
    HadronPhysicsCRMC_FTFP_BERT & operator=( const HadronPhysicsCRMC_FTFP_BERT &right ) = delete;

  protected:
    virtual void Neutron() override;
    virtual void Proton() override;
    virtual void Pion() override;
    virtual void Kaon() override;

  private:
    G4int fModel;                                            // 0:EPOS-LHC, 1:EPOS-1.99, 2:QGSJET:01, 6:SIBYLL-2.3,
    static const std::array< std::string, 13 > fModelNames;  // 7:QGSJETII-04, 11:QGSJETII-03, 12:DPMJET-3.06
    G4double fMinCRMC;
    G4double fMaxFTFP;
};

#endif
