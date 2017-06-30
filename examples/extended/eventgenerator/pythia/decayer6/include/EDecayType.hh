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
// $Id: EDecayType.hh 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/EDecayType.hh
/// \brief Definition of the EDecayType enumeration

#ifndef E_DECAY_TYPE_H
#define E_DECAY_TYPE_H 

/// Enum of decay mode types
///
/// According to EDecayType enum in TPythia6Decayer class in Root:
/// http://root.cern.ch/
/// see http://root.cern.ch/root/License.html
/// ----------------------------------------------------------------------------

enum EDecayType
{
   kSemiElectronic,    
   kDiElectron,
   kSemiMuonic,
   kDiMuon,
   kBJpsiDiMuon,
   kBJpsiDiElectron,
   kBPsiPrimeDiMuon,
   kBPsiPrimeDiElectron,
   kPiToMu,
   kKaToMu,
   kNoDecay,
   kHadronicD,
   kOmega,
   kPhiKK,
   kAll,
   kNoDecayHeavy,
   kHardMuons,
   kBJpsi,
   kWToMuon,
   kWToCharm,
   kWToCharmToMuon,
   kZDiMuon,
   kMaxDecay
};

// ----------------------------------------------------------------------------

#endif
