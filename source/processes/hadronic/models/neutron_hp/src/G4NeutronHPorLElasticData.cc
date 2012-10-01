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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//

#include "G4NeutronHPorLElasticData.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4NeutronHPData.hh"

#include "G4PhysicsVector.hh"

G4NeutronHPorLElasticData::G4NeutronHPorLElasticData()
{
   SetMinKinEnergy( 0*MeV );                                   
   SetMaxKinEnergy( 20*MeV );                                   

   ke_cache = 0.0;
   xs_cache = 0.0;
   element_cache = NULL;
   material_cache = NULL;
//   BuildPhysicsTable(*G4Neutron::Neutron());
}

G4NeutronHPorLElasticData::~G4NeutronHPorLElasticData()
{
//  delete theCrossSections;
}

G4bool G4NeutronHPorLElasticData::IsIsoApplicable( const G4DynamicParticle* dp , 
                                                G4int /*Z*/ , G4int /*A*/ ,
                                                const G4Element* element ,
                                                const G4Material* /*mat*/ )
{
   G4double eKin = dp->GetKineticEnergy();
   if ( eKin > GetMaxKinEnergy() 
     || eKin < GetMinKinEnergy() 
     || dp->GetDefinition() != G4Neutron::Neutron() ) return false;                                   
   if ( unavailable_elements->find( element->GetName() ) != unavailable_elements->end() ) return false;

   return true;
}

G4double G4NeutronHPorLElasticData::GetIsoCrossSection( const G4DynamicParticle* dp ,
                                   G4int /*Z*/ , G4int /*A*/ ,
                                   const G4Isotope* /*iso*/  ,
                                   const G4Element* element ,
                                   const G4Material* material )
{
   if ( dp->GetKineticEnergy() == ke_cache && element == element_cache &&  material == material_cache ) return xs_cache;

   ke_cache = dp->GetKineticEnergy();
   element_cache = element;
   material_cache = material;
   G4double xs = GetCrossSection( dp , element , material->GetTemperature() );
   xs_cache = xs;
   return xs;
   //return GetCrossSection( dp , element , material->GetTemperature() );
}


G4NeutronHPorLElasticData::G4NeutronHPorLElasticData( G4NeutronHPChannel* pChannel , std::set< G4String >* pSet )
:G4VCrossSectionDataSet("NeutronHPorLElasticXS")
{
   theElasticChannel = pChannel;
   unavailable_elements = pSet;   

   SetMinKinEnergy( 0*MeV );                                   
   SetMaxKinEnergy( 20*MeV );                                   

   ke_cache = 0.0;
   xs_cache = 0.0;
   element_cache = NULL;
   material_cache = NULL;
}

 
/*
G4bool G4NeutronHPorLElasticData::IsApplicable(const G4DynamicParticle*aP, const G4Element* anElement)
{
   G4bool result = true;
   G4double eKin = aP->GetKineticEnergy();
   if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
   if ( unavailable_elements->find( anElement->GetName() ) != unavailable_elements->end() ) result = false;
   return result;
}
*/

void G4NeutronHPorLElasticData::BuildPhysicsTable( const G4ParticleDefinition& aP )
{
   if( &aP!=G4Neutron::Neutron() ) 
      throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
}
 


void G4NeutronHPorLElasticData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
//  G4cout << "G4NeutronHPorLElasticData::DumpPhysicsTable still to be implemented"<<G4endl;
}



#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"

G4double G4NeutronHPorLElasticData::
GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double aT)
{

  //G4cout << "Choice G4NeutronHPorLElasticData for element " << anE->GetName() << G4endl;
  G4double result = 0;
//  G4bool outOfRange;
  G4int index = anE->GetIndex();

  // prepare neutron
  G4double eKinetic = aP->GetKineticEnergy();
  G4ReactionProduct theNeutron( aP->GetDefinition() );
  theNeutron.SetMomentum( aP->GetMomentum() );
  theNeutron.SetKineticEnergy( eKinetic );

  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eps = 0.0001;
  G4double theA = anE->GetN();
  G4double theZ = anE->GetZ();
  G4double eleMass; 
  eleMass = ( G4NucleiProperties::GetNuclearMass(static_cast<G4int>(theA+eps), static_cast<G4int>(theZ+eps))
	     ) / G4Neutron::Neutron()->GetPDGMass();
  
  G4ReactionProduct boosted;
  G4double aXsection;
  
  // MC integration loop
  G4int counter = 0;
  G4double buffer = 0;
  G4int size = G4int(std::max(10., aT/60*kelvin));
  G4ThreeVector neutronVelocity = 1./G4Neutron::Neutron()->GetPDGMass()*theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();
  while(counter == 0 || std::abs(buffer-result/std::max(1,counter)) > 0.03*buffer)
  {
    if(counter) buffer = result/counter;
    while (counter<size)
    {
      counter ++;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus(eleMass, aT);
      boosted.Lorentz(theNeutron, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      //aXsection = (*((*theCrossSections)(index))).GetValue(theEkin, outOfRange);
      aXsection = theElasticChannel[index].GetXsec( theEkin );
      // velocity correction.
      G4ThreeVector targetVelocity = 1./aThermalNuc.GetMass()*aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity-neutronVelocity).mag()/neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  return result;
}
