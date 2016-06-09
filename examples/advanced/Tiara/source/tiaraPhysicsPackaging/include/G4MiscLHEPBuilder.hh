//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4MiscLHEPBuilder_h
#define G4MiscLHEPBuilder_h 1

#include "globals.hh"
#include "G4LElastic.hh"

#include "G4HadronElasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
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

#include "G4LEAntiProtonInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

// High-energy Models

#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"

class G4MiscLHEPBuilder 
{
  public: 
    G4MiscLHEPBuilder();
    virtual ~G4MiscLHEPBuilder();

  public: 
    void Build();

  private:
 
   G4LElastic* theElasticModel;

   // anti-proton
   G4HadronElasticProcess theAntiProtonElasticProcess;
   G4AntiProtonInelasticProcess theAntiProtonInelastic;
   G4LEAntiProtonInelastic* theLEAntiProtonModel;
   G4HEAntiProtonInelastic* theHEAntiProtonModel;
    
   // anti-neutron
   G4HadronElasticProcess theAntiNeutronElasticProcess;
   G4AntiNeutronInelasticProcess  theAntiNeutronInelastic;
   G4LEAntiNeutronInelastic* theLEAntiNeutronModel;
   G4HEAntiNeutronInelastic* theHEAntiNeutronModel;
   
   // Lambda
   G4HadronElasticProcess theLambdaElasticProcess;
   G4LambdaInelasticProcess  theLambdaInelastic;
   G4LELambdaInelastic*  theLELambdaModel;
   G4HELambdaInelastic*  theHELambdaModel;
  
   // AntiLambda
   G4HadronElasticProcess theAntiLambdaElasticProcess;
   G4AntiLambdaInelasticProcess  theAntiLambdaInelastic;
   G4LEAntiLambdaInelastic*  theLEAntiLambdaModel;
   G4HEAntiLambdaInelastic*  theHEAntiLambdaModel;
  
   // SigmaMinus
   G4HadronElasticProcess theSigmaMinusElasticProcess;
   G4SigmaMinusInelasticProcess  theSigmaMinusInelastic;
   G4LESigmaMinusInelastic*  theLESigmaMinusModel;
   G4HESigmaMinusInelastic*  theHESigmaMinusModel;
  
   // AntiSigmaMinus
   G4HadronElasticProcess theAntiSigmaMinusElasticProcess;
   G4AntiSigmaMinusInelasticProcess  theAntiSigmaMinusInelastic;
   G4LEAntiSigmaMinusInelastic*  theLEAntiSigmaMinusModel;
   G4HEAntiSigmaMinusInelastic*  theHEAntiSigmaMinusModel;
   
   // SigmaPlus
   G4HadronElasticProcess theSigmaPlusElasticProcess;
   G4SigmaPlusInelasticProcess  theSigmaPlusInelastic;
   G4LESigmaPlusInelastic*  theLESigmaPlusModel;
   G4HESigmaPlusInelastic*  theHESigmaPlusModel;
  
   // AntiSigmaPlus
   G4HadronElasticProcess theAntiSigmaPlusElasticProcess;
   G4AntiSigmaPlusInelasticProcess  theAntiSigmaPlusInelastic;
   G4LEAntiSigmaPlusInelastic*  theLEAntiSigmaPlusModel;
   G4HEAntiSigmaPlusInelastic*  theHEAntiSigmaPlusModel;
  
   // XiZero
   G4HadronElasticProcess theXiZeroElasticProcess;
   G4XiZeroInelasticProcess  theXiZeroInelastic;
   G4LEXiZeroInelastic*  theLEXiZeroModel;
   G4HEXiZeroInelastic*  theHEXiZeroModel;
  
   // AntiXiZero
   G4HadronElasticProcess theAntiXiZeroElasticProcess;
   G4AntiXiZeroInelasticProcess  theAntiXiZeroInelastic;
   G4LEAntiXiZeroInelastic*  theLEAntiXiZeroModel;
   G4HEAntiXiZeroInelastic*  theHEAntiXiZeroModel;
  
   // XiMinus
   G4HadronElasticProcess theXiMinusElasticProcess;
   G4XiMinusInelasticProcess  theXiMinusInelastic;
   G4LEXiMinusInelastic*  theLEXiMinusModel;
   G4HEXiMinusInelastic*  theHEXiMinusModel;

   // AntiXiMinus
   G4HadronElasticProcess theAntiXiMinusElasticProcess;
   G4AntiXiMinusInelasticProcess  theAntiXiMinusInelastic;
   G4LEAntiXiMinusInelastic*  theLEAntiXiMinusModel;
   G4HEAntiXiMinusInelastic*  theHEAntiXiMinusModel;
  
   // OmegaMinus
   G4HadronElasticProcess theOmegaMinusElasticProcess;
   G4OmegaMinusInelasticProcess  theOmegaMinusInelastic;
   G4LEOmegaMinusInelastic*  theLEOmegaMinusModel;
   G4HEOmegaMinusInelastic*  theHEOmegaMinusModel;
   
   // AntiOmegaMinus
   G4HadronElasticProcess theAntiOmegaMinusElasticProcess;
   G4AntiOmegaMinusInelasticProcess  theAntiOmegaMinusInelastic;
   G4LEAntiOmegaMinusInelastic*  theLEAntiOmegaMinusModel;
   G4HEAntiOmegaMinusInelastic*  theHEAntiOmegaMinusModel;
   
};
// 2002 by J.P. Wellisch

#endif
