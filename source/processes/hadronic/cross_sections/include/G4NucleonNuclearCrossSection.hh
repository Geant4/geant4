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
// author: Vladimir.Grichine@cern.ch
//
// Implements data from: Barashenkov V.S., Nucleon-Nucleus Cross Section,
// Preprint JINR P2-89-770, p. 12, Dubna 1989 (scanned version from KEK)
// Based on G. Folger version of G4PiNuclearCrossSection class
//
// Modified:
// 05.03.07 V.Ivanchenko - add IfZAApplicable, remove "debug"
// 06.03.07 V.Ivanchenko - add GetElasticCrossSection for combined dataset
//

#ifndef G4NucleonNuclearCrossSection_h
#define G4NucleonNuclearCrossSection_h

#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"

#include "globals.hh"
#include "G4PiData.hh"
#include "G4HadTmpUtil.hh"

class G4NucleonNuclearCrossSection : public G4VCrossSectionDataSet
{
public:
  
  G4NucleonNuclearCrossSection();
  virtual ~G4NucleonNuclearCrossSection();
    
  static const char* Default_Name() {return "G4NucleonNuclearCrossSection";}

  virtual G4bool IsElementApplicable(const G4DynamicParticle* aParticle,
				     G4int Z, 
				     const G4Material* mat = 0);

  virtual G4double GetElementCrossSection(const G4DynamicParticle* aParticle, 
					  G4int Z, 
					  const G4Material* mat = 0);

  virtual void CrossSectionDescription(std::ostream&) const;

  inline G4double GetElasticCrossSection(const G4DynamicParticle* aParticle, 
					 G4int Z);

  inline G4double GetTotalXsc()  { return fTotalXsc;   };
  inline G4double GetElasticXsc(){ return fElasticXsc; };
  
private:

  G4double Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2);

// add Hydrogen from PDG group.

static const G4double e1[44];

static const G4double he_m_t[44];
static const G4double he_m_in[44];
static const G4double he_p_in[44];

static const G4double be_m_t[44];
static const G4double be_m_in[44];
static const G4double be_p_in[44];

static const G4double c_m_t[44];
static const G4double c_m_in[44];
static const G4double c_p_in[44];


static const G4double e2[44];

static const G4double n_m_t[44];
static const G4double n_m_in[44];
static const G4double n_p_in[44];

static const G4double o_m_t[44];
static const G4double o_m_in[44];
static const G4double o_p_in[44];

static const G4double na_m_t[44];
static const G4double na_m_in[44];
static const G4double na_p_in[44];


static const G4double e3[45];

// static const G4double e3_1[31];

static const G4double al_m_t[45];
static const G4double al_m_in[45];
static const G4double al_p_in[45];

static const G4double si_m_t[45];
static const G4double si_m_in[45];
static const G4double si_p_in[45];

static const G4double ca_m_t[45];
static const G4double ca_m_in[45];
static const G4double ca_p_in[45];


static const G4double e4[47];

static const G4double fe_m_t[47];
static const G4double fe_m_in[47];
static const G4double fe_p_in[47];

static const G4double cu_m_t[47];
static const G4double cu_m_in[47];
static const G4double cu_p_in[47];

static const G4double mo_m_t[47];
static const G4double mo_m_in[47];
static const G4double mo_p_in[47];


static const G4double e5[48];

static const G4double cd_m_t[48];
static const G4double cd_m_in[48];
static const G4double cd_p_in[48];

static const G4double sn_m_t[48];
static const G4double sn_m_in[48];
static const G4double sn_p_in[48];

static const G4double w_m_t[48];
static const G4double w_m_in[48];
static const G4double w_p_in[48];

static const G4double e6[46];

// static const G4double e7[46];

static const G4double pb_m_t[46];
static const G4double pb_m_in[46];
static const G4double pb_p_in[46];

static const G4double u_m_t[46];
static const G4double u_m_in[46];
static const G4double u_p_in[46];

// vectors for treatment

std::vector< G4int >     theZ;
std::vector< G4PiData* > thePipData;
std::vector< G4PiData* > thePimData;

  // cross sections

  G4double fTotalXsc;
  G4double fElasticXsc;

  // particles
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theNeutron;

};

inline
G4double G4NucleonNuclearCrossSection::GetElasticCrossSection(
         const G4DynamicParticle* dp, G4int Z)
{
  GetElementCrossSection(dp, Z);
  return fElasticXsc;
}

#endif
