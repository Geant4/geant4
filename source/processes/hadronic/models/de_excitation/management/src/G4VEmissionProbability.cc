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
//
// $Id: G4VEmissionProbability.cc 66241 2012-12-13 18:34:42Z gunter $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"

G4VEmissionProbability::G4VEmissionProbability()
  :OPTxs(3),useSICB(false),LevelDensity(0.1) 
{
  fG4pow = G4Pow::GetInstance();
  fPairCorr = G4PairingCorrection::GetInstance();
}

G4VEmissionProbability::~G4VEmissionProbability() 
{}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  OPTxs = param->GetDeexModelType();
  LevelDensity = param->GetLevelDensity();
}
