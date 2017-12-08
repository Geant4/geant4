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
// $Id: G4HadronPhysicsINCLXX.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsINCLXX
//
// Author: 2011 P. Kaitaniemi
//
// Modified:
// 19.07.2017 A. Dotti: Refactor code, following FTFP_BERT
// 22.05.2014 D. Mancusi: Extend INCL++ to 20 GeV
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
// 01.03.2013 D. Mancusi: Rename to G4HadronPhysicsINCLXX and introduce
//                        parameters for FTFP and NeutronHP
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 27.11.2011 P.Kaitaniemi: Created using QGSP_INCL_ABLA as a template
//
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsINCLXX_h
#define G4HadronPhysicsINCLXX_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4Cache.hh"
/**
 * Build hadronic physics using INCL++, high-energy models (QGSP or FTFP) and
 * possibly NeutronHP.
 *
 * @see G4INCLXXProtonBuilder
 * @see G4INCLXXNeutronBuilder
 * @see G4INCLXXPionBuilder
 * @see G4IonINCLXXBuilder
 */

class G4VCrossSectionDataSet;
class G4ComponentGGHadronNucleusXsc;


class G4HadronPhysicsINCLXX : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsINCLXX(G4int verbose =1);
    G4HadronPhysicsINCLXX(const G4String& name, const G4bool quasiElastic=true, const G4bool neutronHP=false, const G4bool ftfp=false);
    virtual ~G4HadronPhysicsINCLXX();

  public: 
    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
    virtual void TerminateWorker() override;

    void SetQuasiElastic(G4bool value) {QuasiElastic = value;}; 

  protected:
    virtual void CreateModels();
    virtual void Neutron();
    virtual void Proton();
    virtual void Pion();
    virtual void Kaon();
    virtual void Others();
    //This contains extra configurataion specific to this PL
    virtual void ExtraConfiguration();

    G4bool QuasiElastic;
    G4bool withNeutronHP;
    G4bool withFTFP;
    //Thread-private data write them here to delete them
    G4VectorCache<G4VCrossSectionDataSet*> xs_ds;
    G4Cache<G4ComponentGGHadronNucleusXsc*> xs_k;
};

#endif

