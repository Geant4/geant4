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
// $Id: G4DNAElectronElasticScatteringInWater.hh,v 1.6 2006/10/14 11:03:03 pia Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

// Nucl. Instr. Meth. 155 (1978) 145-156
// J. Phys. D: Appl. Phys. 33 (2000) 932-944
// Phys. Med. Biol. 45 (2000) 3171-3194

#ifndef   G4DNAELECTRONELASTICSCATTERINGINWATER_HH
#define  G4DNAELECTRONELASTICSCATTERINGINWATER_HH 1
 
#include "G4VDNAProcessInWater.hh"
 
// TotalCrossSectionPolicy must provide:
//  - [protected] const G4ParticleDefinition * IncomingParticleDefinition(void)
//  - [protected] G4double TotalCrossSection(G4double k, G4int z)
//  - [protected] void BuildTotalCrossSection(void)
 
// FinalStatesPolicy must provide:
//  - [protected] G4double RandomizeCosTheta(G4double k, G4int z)
//  - [protected] G4bool KillIncomingParticle(G4double k)
//  - [protected] void BuildFinalStatesData(void)
 
template<typename TotalCrossSectionPolicy, typename FinalStatesPolicy>
class G4DNAElectronElasticScatteringInWater : public G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>
{
public:

  G4DNAElectronElasticScatteringInWater(const G4String & name): G4VDNAProcessInWater<TotalCrossSectionPolicy, FinalStatesPolicy>(name) {}

  virtual ~G4DNAElectronElasticScatteringInWater() {}
 
  virtual G4VParticleChange* PostStepDoIt(const G4Track & aTrack, const G4Step & aStep);

private:
  // Hides default constructor and assignment operator as private 
  G4DNAElectronElasticScatteringInWater(const G4DNAElectronElasticScatteringInWater & copy);
  G4DNAElectronElasticScatteringInWater & operator=(const G4DNAElectronElasticScatteringInWater & right);
};
 
#include "G4DNAElectronElasticScatteringInWater.icc"
#endif /* G4DNAELECTRONELASTICSCATTERINGINWATER_HH */

