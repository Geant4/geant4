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
// $Id: G4BGGNucleonInelasticXS.hh,v 1.2 2008/12/01 16:50:23 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BGGNucleonInelasticXS
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.03.2007
// Modifications:
//
//
// Class Description:
//
// Wrapper of proton and neutron inelastic cross-sections using Barashenkov 
// parametersation below 100 GeV and Glauber-Gribov model above
//
// -------------------------------------------------------------------
//

#ifndef G4BGGNucleonInelasticXS_h
#define G4BGGNucleonInelasticXS_h

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4GlauberGribovCrossSection;
class G4NucleonNuclearCrossSection;
class G4HadronNucleonXsc;

class G4BGGNucleonInelasticXS : public G4VCrossSectionDataSet
{
public:

  G4BGGNucleonInelasticXS (const G4ParticleDefinition*);

  virtual ~G4BGGNucleonInelasticXS();
   
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

  G4double CoulombFactor(G4double kinEnergy, G4double A);

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4bool G4BGGNucleonInelasticXS::IsApplicable(const G4DynamicParticle*, 
					     const G4Element*)
{
  return true;
  //  return IsZAApplicable(dp, elm->GetZ(), elm->GetN());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4bool G4BGGNucleonInelasticXS::IsZAApplicable(const G4DynamicParticle*, 
					       G4double /*Z*/, G4double/* A*/)
{
  return true;
  //  return (dp->GetDefinition() == particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4BGGNucleonInelasticXS::GetCrossSection(const G4DynamicParticle* dp, 
						  const G4Element* elm, 
						  G4double temp)
{
  return GetIsoZACrossSection(dp, elm->GetZ(), elm->GetN(), temp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
