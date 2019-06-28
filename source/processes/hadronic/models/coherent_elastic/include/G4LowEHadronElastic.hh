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
// Geant4 Header : G4LowEHadronElastic
//
// Author : V.Ivanchenko 10 May 2019
//  
//
// Class Description:
//
// Elastic scattering in resonance energy region 
//

#ifndef G4LowEHadronElastic_h
#define G4LowEHadronElastic_h 1
 
#include "globals.hh"
#include "G4HadronElastic.hh"

class G4LowEHadronElastic : public G4HadronElastic
{
public:

  explicit G4LowEHadronElastic();

  ~G4LowEHadronElastic() override;
 
  G4double SampleInvariantT(const G4ParticleDefinition* p, 
			    G4double plab, G4int Z, G4int A) override;

private:

  G4bool IsResonanseScattering(const G4ParticleDefinition* p, 
			       G4double plab, G4int Z, G4int A);

  G4double plabLowLimit;
  G4double plabHighLimit;
};

#endif
