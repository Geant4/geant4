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
#ifndef ParticleCodeMap_hh
#define ParticleCodeMap_hh

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4SigmaMinus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4SigmaPlus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4XiMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiOmegaMinus.hh"

template <class T> struct ParticleCodeMap;

template<> struct ParticleCodeMap<G4PionPlus>      {  enum {PDG_CODE =   211}; };
template<> struct ParticleCodeMap<G4PionMinus>     {  enum {PDG_CODE =  -211}; };
template<> struct ParticleCodeMap<G4KaonPlus>      {  enum {PDG_CODE =   321}; };
template<> struct ParticleCodeMap<G4KaonMinus>     {  enum {PDG_CODE =  -321}; };
template<> struct ParticleCodeMap<G4Proton>        {  enum {PDG_CODE =  2212}; };
template<> struct ParticleCodeMap<G4AntiProton>    {  enum {PDG_CODE = -2212}; };
template<> struct ParticleCodeMap<G4SigmaMinus>    {  enum {PDG_CODE =  3112}; };
template<> struct ParticleCodeMap<G4AntiSigmaMinus>{  enum {PDG_CODE = -3112}; };
template<> struct ParticleCodeMap<G4SigmaPlus>     {  enum {PDG_CODE =  3222}; };
template<> struct ParticleCodeMap<G4AntiSigmaPlus> {  enum {PDG_CODE = -3222}; };
template<> struct ParticleCodeMap<G4XiMinus>       {  enum {PDG_CODE =  3312}; };
template<> struct ParticleCodeMap<G4AntiXiMinus>   {  enum {PDG_CODE = -3312}; };
template<> struct ParticleCodeMap<G4OmegaMinus>    {  enum {PDG_CODE =  3334}; };
template<> struct ParticleCodeMap<G4AntiOmegaMinus>{  enum {PDG_CODE = -3334}; };

#endif

	  
