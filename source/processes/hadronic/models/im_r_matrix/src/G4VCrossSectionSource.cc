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

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"
#include "G4ios.hh"
#include "G4Pow.hh"
#include "G4VCrossSectionSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4KineticTrack.hh"
#include "G4CrossSectionVector.hh"
#include "G4CrossSectionSourcePtr.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

G4VCrossSectionSource::G4VCrossSectionSource()
{ }


G4VCrossSectionSource::~G4VCrossSectionSource()
{ }


const G4ParticleDefinition * G4VCrossSectionSource::
FindKeyParticle(const G4KineticTrack& trk1,const G4KineticTrack& trk2) const
{
  const G4ParticleDefinition * result;
  
  const G4ParticleDefinition * p1 = trk1.GetDefinition();
  const G4ParticleDefinition * p2 = trk2.GetDefinition();
  
  if( (p1==G4Proton::Proton() && p2==G4Proton::Proton() ) ||
      (p1==G4Neutron::Neutron() && p2==G4Neutron::Neutron()) )
  {
    result = G4Proton::Proton();
  }
  else if( (p1==G4Neutron::Neutron() && p2==G4Proton::Proton()) ||
           (p2==G4Neutron::Neutron() && p1==G4Proton::Proton()) )
  {
    result = G4Neutron::Neutron();
  }
  else
  {
    throw G4HadronicException(__FILE__, __LINE__, "G4VCrossSectionSource: unklnown particles in FindKeyParticle");
  }
  return result;
}

G4bool G4VCrossSectionSource::operator==(const G4VCrossSectionSource &right) const
{
  return (this == (G4VCrossSectionSource *) &right);
}


G4bool G4VCrossSectionSource::operator!=(const G4VCrossSectionSource &right) const
{
  return (this != (G4VCrossSectionSource *) &right);
}


void G4VCrossSectionSource::Print() const
{
  std::size_t nComponents = 0;
  const G4CrossSectionVector* components = GetComponents();
  if (components)
    {
      nComponents = components->size();
    }
  G4cout << "---- " << this->Name() << " ---- has " << nComponents << " components" <<G4endl;
  for (std::size_t i=0; i<nComponents; ++i)
    {
      G4cout << "-" <<  this->Name() << " - Component " << i << ": " <<G4endl;

      G4CrossSectionSourcePtr componentPtr = (*components)[i];
      G4VCrossSectionSource* component = componentPtr();
      component->Print();
    }
}


void G4VCrossSectionSource::PrintAll(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
  G4double sigma = CrossSection(trk1,trk2) / millibarn;
  G4cout << "---- " << Name() << ": "
	 << "Ecm = " << sqrtS / GeV << " GeV -  " 
	 << " Cross section = " << sigma << " mb "
	 << G4endl;

  std::size_t nComponents = 0;
  const G4CrossSectionVector* components = GetComponents();
  if (components != 0)
    {
      nComponents = components->size();
    }
  for (std::size_t i=0; i<nComponents; ++i)
    {
      G4cout << "* Component " << i << ": ";
      G4CrossSectionSourcePtr componentPtr = (*components)[i];
      G4VCrossSectionSource* component = componentPtr();
      component->PrintAll(trk1,trk2);
    }
}


G4bool G4VCrossSectionSource::InLimits(G4double e, G4double eLow, G4double eHigh) const
{
  G4bool answer = false;
  if (e >= eLow && e <= eHigh) answer = true;
  return answer;
}

G4double G4VCrossSectionSource::LowLimit() const
{
  return 0.; 
}


G4double G4VCrossSectionSource::HighLimit() const
{
  return DBL_MAX; 
}

G4bool G4VCrossSectionSource::IsValid(G4double e) const
{
  G4bool answer = false;
  if (e >= LowLimit() && e <= HighLimit()) answer = true;
  return answer;
}

const G4ParticleDefinition* G4VCrossSectionSource::FindLightParticle(const G4KineticTrack& trk1, 
								     const G4KineticTrack& trk2) const
{
  G4double mass1 = trk1.GetDefinition()->GetPDGMass();
  G4double mass2 = trk2.GetDefinition()->GetPDGMass();
  if (mass1 < mass2)
    {
      return trk1.GetDefinition();
    }
  else
    {
      return trk2.GetDefinition();
    }
}


G4double G4VCrossSectionSource::FcrossX(G4double e, G4double e0, 
                     G4double sigma, G4double eParam, G4double power) const
{
  G4double result = 0.;

  G4double denom = eParam*eParam + (e-e0)*(e-e0);
  if (denom > 0.) 
  {
    G4double value = (2.* eParam * sigma * (e-e0) / denom) * G4Pow::GetInstance()->powA(((e0 + eParam) / e), power);
    result = std::max(0., value);
  }
  return result;
}     




