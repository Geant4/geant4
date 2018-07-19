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
// GEANT4 tag $Name: $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HyperonFTFPBuilder
//
// Author: 2012 G.Folger
//
// Modified:
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
#ifndef G4HyperonFTFPBuilder_h
#define G4HyperonFTFPBuilder_h 1

#include "G4PhysicsBuilderInterface.hh"
#include "globals.hh"

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

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4ChipsHyperonInelasticXS.hh"


class G4HyperonFTFPBuilder : public G4PhysicsBuilderInterface
{
  public: 
    G4HyperonFTFPBuilder();
    virtual ~G4HyperonFTFPBuilder();

    virtual void Build() final override;

  private:
 
    G4TheoFSGenerator * HyperonFTFP;
    G4TheoFSGenerator * AntiHyperonFTFP;
    G4GeneratorPrecompoundInterface * theCascade;
    G4FTFModel * theStringModel;
    G4ExcitedStringDecay * theStringDecay;
    G4LundStringFragmentation * theLund;
    G4CascadeInterface * theBertini;
    G4LambdaInelasticProcess*  theLambdaInelastic;
    G4AntiLambdaInelasticProcess*  theAntiLambdaInelastic;
    G4SigmaMinusInelasticProcess*  theSigmaMinusInelastic;
    G4AntiSigmaMinusInelasticProcess*  theAntiSigmaMinusInelastic;
    G4SigmaPlusInelasticProcess*  theSigmaPlusInelastic;
    G4AntiSigmaPlusInelasticProcess*  theAntiSigmaPlusInelastic;
    G4XiZeroInelasticProcess*  theXiZeroInelastic;
    G4AntiXiZeroInelasticProcess*  theAntiXiZeroInelastic;
    G4XiMinusInelasticProcess*  theXiMinusInelastic;
    G4AntiXiMinusInelasticProcess*  theAntiXiMinusInelastic;
    G4OmegaMinusInelasticProcess*  theOmegaMinusInelastic;
  G4AntiOmegaMinusInelasticProcess*  theAntiOmegaMinusInelastic;
  
  //  G4QHadronInelasticDataSet * theCHIPSInelastic;
  G4VCrossSectionDataSet* theCHIPSInelastic;
    G4bool wasActivated;
};
#endif
