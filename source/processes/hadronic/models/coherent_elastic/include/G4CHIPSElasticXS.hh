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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4CHIPSElasticXS
//
// Author  Ivantchenko, Geant4, 23-SEP-2010
//
// Modifications:
// 31-05-2011 V.Uzhinsky added anti-baryons, Pi+, Pi-, K+, K- cross sections
// 23-08-2011 V.Ivanchenko migration to new design and cleanup
//

// Class Description:
// This is a base class for nucleon elastic cross section based on
// CHIPS parameterisation, the class is extracted from 
// G4UHadronElasticModel
// Class Description - End
 
#ifndef G4CHIPSElasticXS_h
#define G4CHIPSElasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4VQCrossSection;

class G4CHIPSElasticXS : public G4VCrossSectionDataSet
{
public: 

  G4CHIPSElasticXS();

  virtual ~G4CHIPSElasticXS();

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, 
			 G4int /*Z*/, G4int /*A*/,
			 const G4Element*, const G4Material*);

  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, 
			      G4int /*Z*/, G4int /*A*/, 
			      const G4Isotope*, const G4Element*, 
			      const G4Material*);

  virtual void Description() const;

private: 

  G4CHIPSElasticXS & operator=(const G4CHIPSElasticXS &right);
  G4CHIPSElasticXS(const G4CHIPSElasticXS&);
  
  G4VQCrossSection* pCManager;
  G4VQCrossSection* nCManager;

  G4VQCrossSection* PBARxsManager;   
  G4VQCrossSection* PIPxsManager;    
  G4VQCrossSection* PIMxsManager;    
  G4VQCrossSection* KPxsManager;     
  G4VQCrossSection* KMxsManager;     

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;
  const G4ParticleDefinition* theParticle;

  G4double thEnergy;
  G4double lowestEnergy;

  G4int  pPDG;
  G4bool isInitialized;
};

#endif
