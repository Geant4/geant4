//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#ifndef G4PiNuclearCrossSection_h
#define G4PiNuclearCrossSection_h

#include "G4VCrossSectionDataSet.hh"

#include "globals.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4PiData.hh"

class G4PiNuclearCrossSection : public G4VCrossSectionDataSet
{
  public:
  
  G4PiNuclearCrossSection();
  virtual ~G4PiNuclearCrossSection();
  G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Element* anElement)
  {
    G4bool result = false;
    if(aParticle->GetDefinition() == G4PionMinus::PionMinus()) result=true;
    if(aParticle->GetDefinition() == G4PionPlus::PionPlus())   result=true;
    if(anElement->GetZ() == 1) result = false;
    return result;
  }
  G4double GetCrossSection(const G4DynamicParticle* aParticle, 
                           const G4Element* anElement,
                           G4double T=0.);
  void BuildPhysicsTable(const G4ParticleDefinition&) {}
  void DumpPhysicsTable(const G4ParticleDefinition&) {}
  
  private:
  G4double Interpolate(G4int Z1, G4int Z2, G4int Z, G4double x1, G4double x2);

// add Hydrogen from PDG group.

static const G4double e1[38];
static const G4double he_t[38];
static const G4double he_in[38];
static const G4double be_m_t[38];
static const G4double be_m_in[38];
static const G4double be_p_t[24];
static const G4double be_p_in[24];
static const G4double e2[39];
static const G4double c_m_t[39];
static const G4double c_m_in[39];
static const G4double c_p_t[24];
static const G4double c_p_in[24];
static const G4double n_m_t[39];
static const G4double n_m_in[39];
static const G4double n_p_t[27];
static const G4double n_p_in[27];
static const G4double e3[31];
static const G4double o_m_t[31];
static const G4double o_m_in[31];
static const G4double o_p_t[20];
static const G4double o_p_in[20];
static const G4double na_m_t[31];
static const G4double na_m_in[31];
static const G4double na_p_t[22];
static const G4double na_p_in[22];
static const G4double e3_1[31];
static const G4double al_m_t[31];
static const G4double al_m_in[31];
static const G4double al_p_t[21];
static const G4double al_p_in[21];
static const G4double ca_m_t[31];
static const G4double ca_m_in[31];
static const G4double ca_p_t[23];
static const G4double ca_p_in[23];

static const G4double e4[32];
static const G4double fe_m_t[32];
static const G4double fe_m_in[32];
static const G4double fe_p_t[25];
static const G4double fe_p_in[25];
static const G4double cu_m_t[32];
static const G4double cu_m_in[32];
static const G4double cu_p_t[25];
static const G4double cu_p_in[25];
static const G4double e5[34];
static const G4double mo_m_t[34];
static const G4double mo_m_in[34];
static const G4double mo_p_t[27];
static const G4double mo_p_in[27];
static const G4double cd_m_t[34];
static const G4double cd_m_in[34];
static const G4double cd_p_t[28];
static const G4double cd_p_in[28];
static const G4double e6[35];
static const G4double sn_m_t[35];
static const G4double sn_m_in[35];
static const G4double sn_p_t[29];
static const G4double sn_p_in[29];
static const G4double w_m_t[35];
static const G4double w_m_in[35];
static const G4double w_p_t[30];
static const G4double w_p_in[30];
static const G4double e7[35];
static const G4double pb_m_t[35];
static const G4double pb_m_in[35];
static const G4double pb_p_t[30];
static const G4double pb_p_in[30];
static const G4double u_m_t[35];
static const G4double u_m_in[35];
static const G4double u_p_t[30];
static const G4double u_p_in[30];

G4std::vector<G4int> theZ;
G4std::vector<G4PiData *> thePipData;
G4std::vector<G4PiData *> thePimData;

};

#endif
