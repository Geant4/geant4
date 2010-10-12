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
// $Id: G4BGGNucleonElasticXS.hh,v 1.6 2010-10-12 06:02:09 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BGGNucleonElasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// Class Description:
//
// Wrapper of proton and neutron elastic cross-sections using Barashenkov 
// parametersation below 100 GeV and Glauber-Gribov model above
//
// -------------------------------------------------------------------
//

#ifndef G4BGGNucleonElasticXS_h
#define G4BGGNucleonElasticXS_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4HadTmpUtil.hh"


class G4GlauberGribovCrossSection;
class G4NucleonNuclearCrossSection;
class G4HadronNucleonXsc;

class G4BGGNucleonElasticXS : public G4VCrossSectionDataSet
{
public:

  G4BGGNucleonElasticXS (const G4ParticleDefinition*);

  virtual ~G4BGGNucleonElasticXS();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature = 0.);

  virtual
  G4double GetZandACrossSection(const G4DynamicParticle*, G4int /*Z*/,
                                G4int /*A*/, G4double aTemperature = 0.);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&); 

private:

  void Initialise();

  G4double CoulombFactor(G4double kinEnergy, G4int A);

  G4double fGlauberEnergy;  
  G4double fLowEnergy;  
  G4double theGlauberFac[93];
  G4double theCoulombFac[93];

  const G4ParticleDefinition*     particle;
  G4GlauberGribovCrossSection*    fGlauber;
  G4NucleonNuclearCrossSection*   fNucleon;
  G4HadronNucleonXsc*             fHadron;
  G4bool                          isProton;
  G4bool                          isInitialized;
};


inline
G4bool G4BGGNucleonElasticXS::IsApplicable(const G4DynamicParticle*, 
					   const G4Element*)
{
  return true;
}


inline
G4bool G4BGGNucleonElasticXS::IsIsoApplicable(const G4DynamicParticle*,
                                              G4int /*Z*/, G4int/* A*/)
{
  return true;
}


inline
G4double G4BGGNucleonElasticXS::GetCrossSection(const G4DynamicParticle* dp, 
						const G4Element* elm, 
						G4double temp)
{
//  return GetIsoZACrossSection(dp, elm->GetZ(), elm->GetN(), temp);
  G4int Z = G4lrint(elm->GetZ());
  G4int N = G4lrint(elm->GetN());
  return GetZandACrossSection(dp, Z, N, temp);
}

#endif
