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
// $Id: G4ChipsElasticModel.hh 90228 2015-05-21 08:49:57Z gcosmo $
//
// Geant4 Header : G4ChipsElasticModel
//
// Author : V.Ivanchenko 29 June 2009 (redesign old elastic model)
//  
// Modified:
//
// Class Description
// Default model for elastic scattering; GHEISHA algorithm is used 
// Class Description - End

#ifndef G4ChipsElasticModel_h
#define G4ChipsElasticModel_h 1
 
#include "G4HadronElastic.hh"
#include "globals.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ChipsAntiBaryonElasticXS.hh"
#include "G4ChipsPionPlusElasticXS.hh"
#include "G4ChipsPionMinusElasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonMinusElasticXS.hh"


class G4ChipsElasticModel : public G4HadronElastic
{
public:

  G4ChipsElasticModel();

  virtual ~G4ChipsElasticModel();
 
  virtual G4double SampleInvariantT(const G4ParticleDefinition* p, 
				    G4double plab,
				    G4int Z, G4int A);
  virtual void ModelDescription(std::ostream&) const;

private:

  G4ChipsProtonElasticXS* pxsManager;
  G4ChipsNeutronElasticXS* nxsManager;

  G4ChipsAntiBaryonElasticXS* PBARxsManager;
  G4ChipsPionPlusElasticXS* PIPxsManager;
  G4ChipsPionMinusElasticXS* PIMxsManager;
  G4ChipsKaonPlusElasticXS* KPxsManager;
  G4ChipsKaonMinusElasticXS* KMxsManager;

};

#endif
