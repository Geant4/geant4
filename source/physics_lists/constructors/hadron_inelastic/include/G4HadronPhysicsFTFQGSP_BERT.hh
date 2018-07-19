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
// $Id: $
//
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   October 2017
//
// Hadron physics for the new, experimental physics list FTFQGSP_BERT,
// with QGS fragmentation of strings, instead of the Lund string
// fragmentation. Note that the string excitation is still done with FTF,
// exactly as for FTFP_BERT.
// Given that it is an experimental, and perhaps temporary, new type of
// hadron physics, corresponding builders are not created and everything
// is implemented directly in this class.
//----------------------------------------------------------------------------
//
#ifndef G4HadronPhysicsFTFQGSP_BERT_h
#define G4HadronPhysicsFTFQGSP_BERT_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4NeutronRadCapture.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"

#include "G4HadronCaptureProcess.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"

#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiDeuteronInelasticProcess.hh"
#include "G4AntiTritonInelasticProcess.hh"
#include "G4AntiHe3InelasticProcess.hh"
#include "G4AntiAlphaInelasticProcess.hh"

#include "G4ChipsHyperonInelasticXS.hh"


class G4HadronPhysicsFTFQGSP_BERT : public G4VPhysicsConstructor
{
  public: 
    G4HadronPhysicsFTFQGSP_BERT(G4int verbose =1);
    G4HadronPhysicsFTFQGSP_BERT(const G4String& name, G4bool quasiElastic=false);
    virtual ~G4HadronPhysicsFTFQGSP_BERT();

  public: 
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    void CreateModels();

    G4NeutronRadCapture* theNeutronCaptureModel;

    G4PreCompoundModel* thePreEquilib;
    G4GeneratorPrecompoundInterface* theCascade;
    G4FTFModel* theStringModel;
    G4ExcitedStringDecay* theStringDecay;
    G4QGSMFragmentation* theQGSMFragmentation;
    G4ExcitationHandler* theHandler;

    G4TheoFSGenerator* theModel1;
    G4TheoFSGenerator* theModel2;
    G4TheoFSGenerator* theModel3;
    G4CascadeInterface* theBertini1;
    G4CascadeInterface* theBertini2;

    G4HadronCaptureProcess* theNeutronCaptureProcess;

    G4NeutronInelasticProcess*   theNeutronInelastic;
    G4ProtonInelasticProcess*    theProtonInelastic;
    G4PionMinusInelasticProcess* thePionMinusInelastic;
    G4PionPlusInelasticProcess*  thePionPlusInelastic;
    G4KaonMinusInelasticProcess* theKaonMinusInelastic;
    G4KaonPlusInelasticProcess*  theKaonPlusInelastic;
    G4KaonZeroLInelasticProcess* theKaonZeroLInelastic;
    G4KaonZeroSInelasticProcess* theKaonZeroSInelastic;

    G4LambdaInelasticProcess*         theLambdaInelastic;
    G4AntiLambdaInelasticProcess*     theAntiLambdaInelastic;
    G4SigmaMinusInelasticProcess*     theSigmaMinusInelastic;
    G4AntiSigmaMinusInelasticProcess* theAntiSigmaMinusInelastic;
    G4SigmaPlusInelasticProcess*      theSigmaPlusInelastic;
    G4AntiSigmaPlusInelasticProcess*  theAntiSigmaPlusInelastic;
    G4XiZeroInelasticProcess*         theXiZeroInelastic;
    G4AntiXiZeroInelasticProcess*     theAntiXiZeroInelastic;
    G4XiMinusInelasticProcess*        theXiMinusInelastic;
    G4AntiXiMinusInelasticProcess*    theAntiXiMinusInelastic;
    G4OmegaMinusInelasticProcess*     theOmegaMinusInelastic;
    G4AntiOmegaMinusInelasticProcess* theAntiOmegaMinusInelastic;

    G4AntiProtonInelasticProcess*   theAntiProtonInelastic;
    G4AntiNeutronInelasticProcess*  theAntiNeutronInelastic;
    G4AntiDeuteronInelasticProcess* theAntiDeuteronInelastic;
    G4AntiTritonInelasticProcess*   theAntiTritonInelastic;
    G4AntiHe3InelasticProcess*      theAntiHe3Inelastic;
    G4AntiAlphaInelasticProcess*    theAntiAlphaInelastic;

    G4VCrossSectionDataSet* thePiXS;
    G4VCrossSectionDataSet* theKaonXS;
    G4VCrossSectionDataSet* theChipsHyperonInelasticXS;
    G4VCrossSectionDataSet* theAntiNucleonXS;
    G4VCrossSectionDataSet* theNeutronInelasticXS;
    G4VCrossSectionDataSet* theNeutronCaptureXS;
};

#endif

