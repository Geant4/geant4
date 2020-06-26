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
//---------------------------------------------------------------------------
// ClassName: G4HyperonBuilder
// Author: Alberto Ribon
// Date: May 2020
// Description: builder class that collects specific model builders for
//              hyperons and anti-hyperons
// Modified:
//---------------------------------------------------------------------------

#ifndef G4HyperonBuilder_h
#define G4HyperonBuilder_h 1

#include "G4PhysicsBuilderInterface.hh"
#include "globals.hh"
#include "G4VHyperonBuilder.hh"
#include <vector>

class G4LambdaInelasticProcess;
class G4AntiLambdaInelasticProcess;
class G4SigmaMinusInelasticProcess;
class G4AntiSigmaMinusInelasticProcess;
class G4SigmaPlusInelasticProcess;
class G4AntiSigmaPlusInelasticProcess;
class G4XiMinusInelasticProcess;
class G4AntiXiMinusInelasticProcess;
class G4XiZeroInelasticProcess;
class G4AntiXiZeroInelasticProcess;
class G4OmegaMinusInelasticProcess;
class G4AntiOmegaMinusInelasticProcess;


class G4HyperonBuilder : public G4PhysicsBuilderInterface {
  public: 
    G4HyperonBuilder();
    virtual ~G4HyperonBuilder() {}
    virtual void Build() final override;
    virtual void RegisterMe( G4PhysicsBuilderInterface* aB ) final override;
  private:
    G4LambdaInelasticProcess*         theLambdaInelastic;
    G4AntiLambdaInelasticProcess*     theAntiLambdaInelastic;
    G4SigmaMinusInelasticProcess*     theSigmaMinusInelastic;
    G4AntiSigmaMinusInelasticProcess* theAntiSigmaMinusInelastic;
    G4SigmaPlusInelasticProcess*      theSigmaPlusInelastic;
    G4AntiSigmaPlusInelasticProcess*  theAntiSigmaPlusInelastic;
    G4XiMinusInelasticProcess*        theXiMinusInelastic;
    G4AntiXiMinusInelasticProcess*    theAntiXiMinusInelastic;
    G4XiZeroInelasticProcess*         theXiZeroInelastic;
    G4AntiXiZeroInelasticProcess*     theAntiXiZeroInelastic;
    G4OmegaMinusInelasticProcess*     theOmegaMinusInelastic;
    G4AntiOmegaMinusInelasticProcess* theAntiOmegaMinusInelastic;
    std::vector< G4VHyperonBuilder* > theModelCollections;
};

#endif
