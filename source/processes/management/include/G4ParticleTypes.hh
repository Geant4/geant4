// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleTypes.hh,v 1.4 2000-02-25 10:11:18 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      4-th April 1996, G.Cosmo
//     ----------------- G4ParticleTypes.hh -----------------
// Class Description
//  Files including all existing particle definitions GEANT4 classes 
//
// ****************************************************************
// Created, G.Cosmo, 19 September 1996
// Modified H.Kurashige 8 Mar. 1997
// Modified H.Kurashige 18 June 1997
// Moved from global/management, G.Cosmo, 1 December 1998
// ----------------------------------------------------------------

#ifndef G4ParticleTypes_h
#define G4ParticleTypes_h 1

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// Leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Mesons
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"

#include "G4DMesonPlus.hh"
#include "G4DMesonMinus.hh"
#include "G4DMesonZero.hh"
#include "G4AntiDMesonZero.hh"
#include "G4DsMesonPlus.hh"
#include "G4DsMesonMinus.hh"
#include "G4JPsi.hh"

#include "G4BMesonPlus.hh"
#include "G4BMesonMinus.hh"
#include "G4BMesonZero.hh"
#include "G4AntiBMesonZero.hh"
#include "G4BsMesonZero.hh"
#include "G4AntiBsMesonZero.hh"


// Barions
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"

#include "G4LambdacPlus.hh"
#include "G4SigmacPlusPlus.hh"
#include "G4SigmacPlus.hh"
#include "G4SigmacZero.hh"
#include "G4XicPlus.hh"
#include "G4XicZero.hh"
#include "G4OmegacZero.hh"

#include "G4AntiLambdacPlus.hh"
#include "G4AntiSigmacPlusPlus.hh"
#include "G4AntiSigmacPlus.hh"
#include "G4AntiSigmacZero.hh"
#include "G4AntiXicPlus.hh"
#include "G4AntiXicZero.hh"
#include "G4AntiOmegacZero.hh"

// Nuclei
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4He3.hh"
#include "G4Triton.hh"

//ions
#include "G4GenericIon.hh"
#endif

