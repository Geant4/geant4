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
// Split barion (antibarion) into quark and diquark (antidiquark and antiqaurk ) 
// based on prototype, needs clean up of interfaces HPW Feb 1999  
// Numbers verified and errors corrected, HPW Dec 1999

#include "G4BaryonSplitter.hh"
#include "G4ParticleTable.hh"

G4BaryonSplitter::
G4BaryonSplitter()
{
  theBaryons.insert(new G4SPBaryon(G4Proton::Proton()));
  theBaryons.insert(new G4SPBaryon(G4Neutron::Neutron()));
  theBaryons.insert(new G4SPBaryon(G4AntiProton::AntiProton()));
  theBaryons.insert(new G4SPBaryon(G4AntiNeutron::AntiNeutron()));
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(2224))); // Delta++
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(2214))); // Delta+
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(2114))); // Delta0
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(1114))); // Delta-
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(-2224))); // anti Delta++
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(-2214))); // anti Delta+
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(-2114))); // anti Delta0
  theBaryons.insert(new G4SPBaryon(G4ParticleTable::GetParticleTable()->FindParticle(-1114))); // anti Delta-
  theBaryons.insert(new G4SPBaryon(G4Lambda::Lambda()));
  theBaryons.insert(new G4SPBaryon(G4AntiLambda::AntiLambda()));
  theBaryons.insert(new G4SPBaryon(G4SigmaPlus::SigmaPlus()));
  theBaryons.insert(new G4SPBaryon(G4SigmaZero::SigmaZero()));
  theBaryons.insert(new G4SPBaryon(G4SigmaMinus::SigmaMinus()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmaPlus::AntiSigmaPlus()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmaZero::AntiSigmaZero()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmaMinus::AntiSigmaMinus()));
  theBaryons.insert(new G4SPBaryon(G4XiMinus::XiMinus()));
  theBaryons.insert(new G4SPBaryon(G4XiZero::XiZero()));
  theBaryons.insert(new G4SPBaryon(G4AntiXiMinus::AntiXiMinus()));
  theBaryons.insert(new G4SPBaryon(G4AntiXiZero::AntiXiZero()));
  theBaryons.insert(new G4SPBaryon(G4OmegaMinus::OmegaMinus()));
  theBaryons.insert(new G4SPBaryon(G4AntiOmegaMinus::AntiOmegaMinus()));
  theBaryons.insert(new G4SPBaryon(G4LambdacPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiLambdacPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmacPlusPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmacPlusPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmacPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmacPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmacZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmacZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4XicPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiXicPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4XicZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiXicZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4OmegacZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiOmegacZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4Lambdab::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiLambdab::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmabPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmabPlus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmabZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmabZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4SigmabMinus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiSigmabMinus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4XibZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiXibZero::Definition()));
  theBaryons.insert(new G4SPBaryon(G4XibMinus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiXibMinus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4OmegabMinus::Definition()));
  theBaryons.insert(new G4SPBaryon(G4AntiOmegabMinus::Definition()));
}

G4bool G4BaryonSplitter::
SplitBarion(G4int PDGCode, G4int* q_or_qqbar, G4int* qbar_or_qq)
{
  const G4SPBaryon * aBaryon = theBaryons.GetBaryon(G4ParticleTable::GetParticleTable()->FindParticle(PDGCode));

  if(aBaryon==NULL)
  {
    return FALSE;
   } else {
    aBaryon->SampleQuarkAndDiquark(*q_or_qqbar, *qbar_or_qq);
    return TRUE;
  }
}


// Get the splittable baryon for a PDG code.
const G4SPBaryon & G4BaryonSplitter::
GetSPBaryon(G4int PDGCode)
{
  return *theBaryons.GetBaryon(G4ParticleTable::GetParticleTable()->FindParticle(PDGCode));
}


// Find rest diquark in given barion after quark - antiquark annihilation  
G4bool G4BaryonSplitter::
FindDiquark(G4int PDGCode, G4int Quark, G4int* Diquark)
{
  const G4SPBaryon * aBaryon = theBaryons.GetBaryon(G4ParticleTable::GetParticleTable()->FindParticle(PDGCode));
  if(aBaryon)
  {
    aBaryon->FindDiquark(Quark, *Diquark);
    return true;
  }
  return false;
}

