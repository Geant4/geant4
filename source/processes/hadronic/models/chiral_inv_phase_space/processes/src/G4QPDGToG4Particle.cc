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
// $Id: G4QPDGToG4Particle.cc,v 1.1 2009-11-17 10:36:55 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------- G4QG4ToG4Particle singletone class ------------------
//                 by Mikhail Kossov, December 2003.
// PDGCode->G4Particle convertor of the CHIPS Simulation Branch in GEANT4
// ----------------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// Short description: This is a helper class, which converts the PDG-defined
// G4QHadrons of the CHIPS model to the G4 particles, defined by the singetones.
// -----------------------------------------------------------------------------


//#define pdebug

#include "G4QPDGToG4Particle.hh"

G4QPDGToG4Particle::G4QPDGToG4Particle()
{
}

G4QPDGToG4Particle::~G4QPDGToG4Particle() // The map is distructed only in the EndOfJob
{
}

// Returns Pointer to the G4QPDGToG4Particle
G4QPDGToG4Particle* G4QPDGToG4Particle::Get()
//                  =========================
{
  static G4QPDGToG4Particle theMap;        // *** Static body of the G4QPDGToG4Particle ***
  return &theMap;
}

G4ParticleDefinition* G4QPDGToG4Particle::GetParticleDefinition(G4int PDG)
//                    ====================================================
{
  if(!PDG) return 0;
  else if(PDG>0)     // Positive PDG Code
  {
    if(PDG<100)
    {
      if(PDG==22) return G4Gamma::Gamma();
      else if(PDG>10 && PDG<17)
      {
        if(PDG<13)
        {
          if(PDG==11) return G4Electron::Electron();
          else return G4NeutrinoE::NeutrinoE();
        }
        else
        {
          if(PDG<15)
          {
            if(PDG==13) return G4MuonMinus::MuonMinus();
            else return G4NeutrinoMu::NeutrinoMu();
          }
          else
          {
            if(PDG==15) return G4TauMinus::TauMinus();
            else return G4NeutrinoTau::NeutrinoTau();
          }
        }
      }
      else return 0; // @@ Warning can be added
    } // End of the Lepton definition
    else if(PDG<1000)
    {
      if(PDG<420)
      {
        if(PDG<320)
        {
          if(PDG==211) return G4PionPlus::PionPlus();
          else if(PDG==111) return G4PionZero::PionZero();
          else if(PDG==130) return G4KaonZeroLong::KaonZeroLong();
          else if(PDG==221) return G4Eta::Eta();
          else if(PDG==311) return G4KaonZero::KaonZero();
          else return 0; // @@ Warning can be added
        }
        else
        {
          if(PDG==321) return G4KaonPlus::KaonPlus();
          else if(PDG==331) return G4EtaPrime::EtaPrime();
          else if(PDG==310) return G4KaonZeroShort::KaonZeroShort();
          else if(PDG==411) return G4DMesonPlus::DMesonPlus();
          else return 0; // @@ Warning can be added
        }
      }
      else
      {
        if(PDG<500)
        {
          if(PDG==421) return G4DMesonZero::DMesonZero();
          else if(PDG==431) return G4DsMesonPlus::DsMesonPlus();
          else if(PDG==443) return G4JPsi::JPsi();
          else return 0; // @@ Warning can be added
        }
        else
        {
          if(PDG==521) return G4BMesonPlus::BMesonPlus();
          else if(PDG==511) return G4BMesonZero::BMesonZero();
          else if(PDG==531) return G4BsMesonZero::BsMesonZero();
          else return 0; // @@ Warning can be added
        }
      }
    } // Emd of the Meson definition
    else
    {
      if(PDG<3333)
      {
        if(PDG<3211)
        {
          if(PDG<3111)
          {
            if(PDG==2112) return G4Neutron::Neutron();
            else if(PDG==2212) return G4Proton::Proton();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(PDG==3112) return G4SigmaMinus::SigmaMinus();
            else if(PDG==3122) return G4Lambda::Lambda();
            else return 0; // @@ Warning can be added
          }
        }
        else
        {
          if(PDG<3311)
          {
            if(PDG==3222) return G4SigmaPlus::SigmaPlus();
            else if(PDG==3212) return G4SigmaZero::SigmaZero();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(PDG==3312) return G4XiMinus::XiMinus();
            else if(PDG==3322) return G4XiZero::XiZero();
            else return 0; // @@ Warning can be added
          }
        }
      }
      else
      {
        if(PDG<4221)
        {
          if(PDG<4121)
          {
            if(PDG==3334) return G4OmegaMinus::OmegaMinus();
            else if(PDG==4112) return G4SigmacZero::SigmacZero();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(PDG==4122) return G4LambdacPlus::LambdacPlus();
            else if(PDG==4212) return G4SigmacPlus::SigmacPlus();
            else return 0; // @@ Warning can be added
          }
        }
        else
        {
          if(PDG<4231)
          {
            if(PDG==4222) return G4SigmacPlusPlus::SigmacPlusPlus();
            else if(PDG==4232) return G4XicPlus::XicPlus();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(PDG==4132) return G4XicZero::XicZero();
            else if(PDG==4332) return G4OmegacZero::OmegacZero();
            else return 0; // @@ Warning can be added
          }
        }
      }
    } // End of Baryon definition
  } 
  else               // Negative PDG Code
  {
    G4int aPDG=-PDG;
#ifdef pdebug
    G4cout<<"G4QPDGToG4Particle:Antiparticle PDG="<<PDG<<G4endl;
#endif
    if(aPDG<100)
    {
      if(aPDG>10 && aPDG<17)
      {
        if(aPDG<13)
        {
          if(aPDG==11) return G4Positron::Positron();
          else return G4AntiNeutrinoE::AntiNeutrinoE();
        }
        else
        {
          if(aPDG<15)
          {
            if(aPDG==13) return G4MuonPlus::MuonPlus();
            else return G4AntiNeutrinoMu::AntiNeutrinoMu();
          }
          else
          {
            if(aPDG==15) return G4TauPlus::TauPlus();
            else return G4AntiNeutrinoTau::AntiNeutrinoTau();
          }
        }
      }
      else return 0; // @@ Warning can be added
    } // End of the Anti-Lepton definition
    else if(aPDG<1000)
    {
#ifdef pdebug
      G4cout<<"G4QPDGToG4Particle:AntiMesons aPDG="<<aPDG<<G4endl;
#endif
      if(aPDG<420)
      {
#ifdef pdebug
       G4cout<<"G4QPDGToG4Particle:AntiSU(3)Mesons aPDG="<<aPDG<<G4endl;
#endif
        if(aPDG<320)
        {
#ifdef pdebug
          G4cout<<"G4QPDGToG4Particle:AntiPi&KMesons aPDG="<<aPDG<<G4endl;
#endif
          if(aPDG==211) return G4PionMinus::PionMinus();
          else if(aPDG==311) return G4AntiKaonZero::AntiKaonZero();
          else return 0; // @@ Warning can be added
        }
        else
        {
#ifdef pdebug
          G4cout<<"G4QPDGToG4Particle:AntiK&DMesons aPDG="<<aPDG<<G4endl;
#endif
          if(aPDG==321) 
          {
#ifdef pdebug
            G4cout<<"G4QPDGToG4Particle:KaonMinus aPDG="<<aPDG<<G4endl;
#endif
            return G4KaonMinus::KaonMinus();
          }
          else if(aPDG==411) return G4DMesonMinus::DMesonMinus();
          else return 0; // @@ Warning can be added
        }
      }
      else
      {
        if(aPDG<500)
        {
          if(aPDG==421) return G4AntiDMesonZero::AntiDMesonZero();
          else if(aPDG==431) return G4DsMesonMinus::DsMesonMinus();
          else return 0; // @@ Warning can be added
        }
        else
        {
          if(aPDG==521) return G4BMesonMinus::BMesonMinus();
          else if(aPDG==511) return G4AntiBMesonZero::AntiBMesonZero();
          else if(aPDG==531) return G4AntiBsMesonZero::AntiBsMesonZero();
          else return 0; // @@ Warning can be added
        }
      }
    } // Emd of the Anti-Meson definition
    else
    {
      if(aPDG<3333)
      {
        if(aPDG<3211)
        {
          if(aPDG<3111)
          {
            if(aPDG==2112) return G4AntiNeutron::AntiNeutron();
            else if(aPDG==2212) return G4AntiProton::AntiProton();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(aPDG==3112) return G4AntiSigmaMinus::AntiSigmaMinus();
            else if(aPDG==3122) return G4AntiLambda::AntiLambda();
            else return 0; // @@ Warning can be added
          }
        }
        else
        {
          if(aPDG<3311)
          {
            if(aPDG==3222) return G4AntiSigmaPlus::AntiSigmaPlus();
            else if(aPDG==3212) return G4AntiSigmaZero::AntiSigmaZero();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(aPDG==3312) return G4AntiXiMinus::AntiXiMinus();
            else if(aPDG==3322) return G4AntiXiZero::AntiXiZero();
            else return 0; // @@ Warning can be added
          }
        }
      }
      else
      {
        if(aPDG<4221)
        {
          if(aPDG<4121)
          {
            if(aPDG==3334) return G4AntiOmegaMinus::AntiOmegaMinus();
            else if(aPDG==4112) return G4AntiSigmacZero::AntiSigmacZero();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(aPDG==4122) return G4AntiLambdacPlus::AntiLambdacPlus();
            else if(aPDG==4212) return G4AntiSigmacPlus::AntiSigmacPlus();
            else return 0; // @@ Warning can be added
          }
        }
        else
        {
          if(aPDG<4231)
          {
            if(aPDG==4222) return G4AntiSigmacPlusPlus::AntiSigmacPlusPlus();
            else if(aPDG==4232) return G4AntiXicPlus::AntiXicPlus();
            else return 0; // @@ Warning can be added
          }
          else
          {
            if(aPDG==4132) return G4AntiXicZero::AntiXicZero();
            else if(aPDG==4332) return G4AntiOmegacZero::AntiOmegacZero();
            else return 0; // @@ Warning can be added
          }
        }
      }
    } // End of Anti-Baryon definition
  } // End of Anti-particle definition
  return 0;
}

void G4QPDGToG4Particle::DefineAllParticles()
// ==========================================
{
  //======== LEPTONS =========
  G4Gamma::GammaDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  //================ MESONS ===================
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
  G4DMesonPlus::DMesonPlusDefinition();
  G4DMesonMinus::DMesonMinusDefinition();
  G4DMesonZero::DMesonZeroDefinition();
  G4AntiDMesonZero::AntiDMesonZeroDefinition();
  G4DsMesonPlus::DsMesonPlusDefinition();
  G4DsMesonMinus::DsMesonMinusDefinition();
  G4JPsi::JPsiDefinition();
  G4BMesonPlus::BMesonPlusDefinition();
  G4BMesonMinus::BMesonMinusDefinition();
  G4BMesonZero::BMesonZeroDefinition();
  G4AntiBMesonZero::AntiBMesonZeroDefinition();
  G4BsMesonZero::BsMesonZeroDefinition();
  G4AntiBsMesonZero::AntiBsMesonZeroDefinition();
  // ========== BARYONS ==================
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
  G4Lambda::LambdaDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4XiMinus::XiMinusDefinition();
  G4XiZero::XiZeroDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiLambda::AntiLambdaDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4AntiSigmaZero::AntiSigmaZeroDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4AntiXiZero::AntiXiZeroDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();
  G4LambdacPlus::LambdacPlusDefinition();
  G4SigmacPlusPlus::SigmacPlusPlusDefinition();
  G4SigmacPlus::SigmacPlusDefinition();
  G4SigmacZero::SigmacZeroDefinition();
  G4XicPlus::XicPlusDefinition();
  G4XicZero::XicZeroDefinition();
  G4OmegacZero::OmegacZeroDefinition();
  G4AntiLambdacPlus::AntiLambdacPlusDefinition();
  G4AntiSigmacPlusPlus::AntiSigmacPlusPlusDefinition();
  G4AntiSigmacPlus::AntiSigmacPlusDefinition();
  G4AntiSigmacZero::AntiSigmacZeroDefinition();
  G4AntiXicPlus::AntiXicPlusDefinition();
  G4AntiXicZero::AntiXicZeroDefinition();
  G4AntiOmegacZero::AntiOmegacZeroDefinition();
}
