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
// $Id: G4DiffElasticRatio.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4DiffElasticRatio
//
// Author:  V.Grichine 07.11.2014
//
// Modifications:
//
 
//
// Class Description
// This is a class for hadronic cross section ratio
// diffractive/elastic
// Class Description - End

#ifndef G4DiffElasticRatio_h
#define G4DiffElasticRatio_h 1

#include "globals.hh"
#include "G4VCrossSectionRatio.hh"

class G4ParticleDefinition;
class G4ComponentGGHadronNucleusXsc;

class G4DiffElasticRatio: public G4VCrossSectionRatio
{
public: 

  G4DiffElasticRatio(const G4String& nam = "", G4int verb = 0);

  virtual ~G4DiffElasticRatio();

  
  G4double ComputeRatio(const G4ParticleDefinition*,
			G4double kinEnergy, 
			G4int Z, G4int A);

  void SetEnergyThreshold(G4double e){fDDthreshold=e;};
  G4double GetEnergyThreshold(){return fDDthreshold;};

private:

  G4DiffElasticRatio & operator=(const G4DiffElasticRatio &right);
  G4DiffElasticRatio(const G4DiffElasticRatio&);

  G4ComponentGGHadronNucleusXsc* fGGXsc;
  G4double fDDthreshold;
};

#endif
