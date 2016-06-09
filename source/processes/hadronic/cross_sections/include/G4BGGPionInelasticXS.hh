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
// $Id: G4BGGPionInelasticXS.hh,v 1.1 2007/03/13 15:19:30 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BGGPionInelasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 01.10.2003
// Modifications:
//
//
// Class Description:
//
// Wrapper of pi+ and pi- inelastic cross-sections using Barashenkov 
// parametersation below 100 GeV and Glauber-Gribov model above
//
// -------------------------------------------------------------------
//

#ifndef G4BGGPionInelasticXS_h
#define G4BGGPionInelasticXS_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4GlauberGribovCrossSection;
class G4UPiNuclearCrossSection;

class G4BGGPionInelasticXS : public G4VCrossSectionDataSet
{
public:

  G4BGGPionInelasticXS (const G4ParticleDefinition*);

  virtual ~G4BGGPionInelasticXS();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4bool IsZAApplicable(const G4DynamicParticle*, G4double Z, G4double A);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature = 0.);

  virtual
  G4double GetIsoZACrossSection(const G4DynamicParticle*, G4double /*Z*/,
				G4double /*A*/, G4double aTemperature = 0.);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&); 

private:

  void Initialise();

  G4double thEnergy;  // threshold of Glauber model
  G4double theFac[93];

  const G4ParticleDefinition*     particle;
  G4GlauberGribovCrossSection*    fGlauber;
  G4UPiNuclearCrossSection*       fPion;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4bool G4BGGPionInelasticXS::IsApplicable(const G4DynamicParticle* dp, 
					const G4Element*  elm)
{
  return IsZAApplicable(dp, elm->GetZ(), elm->GetN());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4bool G4BGGPionInelasticXS::IsZAApplicable(const G4DynamicParticle* dp, 
					  G4double Z, G4double/* A*/)
{
  return (dp->GetDefinition() == particle && Z > 1.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4BGGPionInelasticXS::GetCrossSection(const G4DynamicParticle* dp, 
					     const G4Element* elm, 
					     G4double temp)
{
  return GetIsoZACrossSection(dp, elm->GetZ(), elm->GetN(), temp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
