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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070618 fix memory leaking by T. Koi
// 071002 enable cross section dump by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
// 081124 Protect invalid read which caused run time errors by T. Koi
// 100729 Add safty for 0 lenght cross sections by T. Ko
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4ParticleHPData.hh"
#include "G4ParticleHPManager.hh"
#include "G4Pow.hh"

G4ParticleHPFissionData::G4ParticleHPFissionData()
:G4VCrossSectionDataSet("NeutronHPFissionXS")
{
   SetMinKinEnergy( 0*MeV );                                   
   SetMaxKinEnergy( 20*MeV );                                   

   theCrossSections = 0;
   onFlightDB = true;
   instanceOfWorker = false;
   if ( G4Threading::IsWorkerThread() ) {
      instanceOfWorker = true;
   }
   element_cache = NULL;
   material_cache = NULL;
   ke_cache = 0.0; 
   xs_cache = 0.0; 
   //BuildPhysicsTable(*G4Neutron::Neutron());
}
   
G4ParticleHPFissionData::~G4ParticleHPFissionData()
{
   if ( theCrossSections != NULL && instanceOfWorker != true ) {
     theCrossSections->clearAndDestroy();
     delete theCrossSections;
     theCrossSections = NULL;
   }
}

G4bool G4ParticleHPFissionData::IsIsoApplicable( const G4DynamicParticle* dp , 
                                                G4int /*Z*/ , G4int /*A*/ ,
                                                const G4Element* /*elm*/ ,
                                                const G4Material* /*mat*/ )
{
   G4double eKin = dp->GetKineticEnergy();
   if ( eKin > GetMaxKinEnergy() 
     || eKin < GetMinKinEnergy() 
     || dp->GetDefinition() != G4Neutron::Neutron() ) return false;                                   

   return true;
}

G4double G4ParticleHPFissionData::GetIsoCrossSection( const G4DynamicParticle* dp ,
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
}

/*
G4bool G4ParticleHPFissionData::IsApplicable(const G4DynamicParticle*aP, const G4Element*)
{
  G4bool result = true;
  G4double eKin = aP->GetKineticEnergy();
  if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
  return result;
}
*/

void G4ParticleHPFissionData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{

   if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() ) {
      G4cout << "Find a flag of \"G4NEUTRONHP_NEGLECT_DOPPLER\"." << G4endl;
      G4cout << "On the fly Doppler broadening will be neglect in the cross section calculation of fission reaction of neutrons (<20MeV)." << G4endl;
      onFlightDB = false;
   } 

  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

   if ( G4Threading::IsWorkerThread() ) {
      theCrossSections = G4ParticleHPManager::GetInstance()->GetFissionCrossSections();
      return;
   }

  size_t numberOfElements = G4Element::GetNumberOfElements();
  //theCrossSections = new G4PhysicsTable( numberOfElements );
   // TKDB
   //if ( theCrossSections == NULL ) theCrossSections = new G4PhysicsTable( numberOfElements );
   if ( theCrossSections == NULL ) 
      theCrossSections = new G4PhysicsTable( numberOfElements );
   else
      theCrossSections->clearAndDestroy();

  // make a PhysicsVector for each element

  static G4ThreadLocal G4ElementTable *theElementTable  = 0 ; if (!theElementTable) theElementTable= G4Element::GetElementTable();
  for( size_t i=0; i<numberOfElements; ++i )
  {
    G4PhysicsVector* physVec = G4ParticleHPData::
      Instance(G4Neutron::Neutron())->MakePhysicsVector((*theElementTable)[i], this);
    theCrossSections->push_back(physVec);
  }

   G4ParticleHPManager::GetInstance()->RegisterFissionCrossSections( theCrossSections );
}

void G4ParticleHPFissionData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

//
// Dump element based cross section
// range 10e-5 eV to 20 MeV
// 10 point per decade
// in barn
//

   G4cout << G4endl;
   G4cout << G4endl;
   G4cout << "Fission Cross Section of Neutron HP"<< G4endl;
   G4cout << "(Pointwise cross-section at 0 Kelvin.)" << G4endl;
   G4cout << G4endl;
   G4cout << "Name of Element" << G4endl;
   G4cout << "Energy[eV]  XS[barn]" << G4endl;
   G4cout << G4endl;

   size_t numberOfElements = G4Element::GetNumberOfElements();
   static G4ThreadLocal G4ElementTable *theElementTable  = 0 ; if (!theElementTable) theElementTable= G4Element::GetElementTable();

   for ( size_t i = 0 ; i < numberOfElements ; ++i )
   {

      G4cout << (*theElementTable)[i]->GetName() << G4endl;

      if ( (*((*theCrossSections)(i))).GetVectorLength() == 0 ) 
      {
         G4cout << "The cross-section data of the fission of this element is not available." << G4endl; 
         G4cout << G4endl; 
         continue;
      }

      G4int ie = 0;

      for ( ie = 0 ; ie < 130 ; ie++ )
      {
         G4double eKinetic = 1.0e-5 * G4Pow::GetInstance()->powA ( 10.0 , ie/10.0 ) *eV;
         G4bool outOfRange = false;

         if ( eKinetic < 20*MeV )
         {
            G4cout << eKinetic/eV << " " << (*((*theCrossSections)(i))).GetValue(eKinetic, outOfRange)/barn << G4endl;
         }

      }

      G4cout << G4endl;
   }

  //G4cout << "G4ParticleHPFissionData::DumpPhysicsTable still to be implemented"<<G4endl;
}

#include "G4NucleiProperties.hh"

G4double G4ParticleHPFissionData::
GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double aT)
{
  G4double result = 0;
  if(anE->GetZ()<88) return result;
  G4bool outOfRange;
  G4int index = anE->GetIndex();

// 100729 TK add safety
if ( ( ( *theCrossSections )( index ) )->GetVectorLength() == 0 ) return result;

  // prepare neutron
  G4double eKinetic = aP->GetKineticEnergy();
  G4ReactionProduct theNeutronRP( aP->GetDefinition() );
  theNeutronRP.SetMomentum( aP->GetMomentum() );
  theNeutronRP.SetKineticEnergy( eKinetic );

  if ( !onFlightDB ) {
     //NEGLECT_DOPPLER
     G4double factor = 1.0;
     if ( eKinetic < aT * k_Boltzmann ) {
        // below 0.1 eV neutrons 
        // Have to do some, but now just igonre.   
        // Will take care after performance check.  
        // factor = factor * targetV;
     }
     return ( (*((*theCrossSections)(index))).GetValue(eKinetic, outOfRange) )* factor; 
  }

  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eps = 0.0001;
  G4double theA = anE->GetN();
  G4double theZ = anE->GetZ();
  G4double eleMass; 
  eleMass = ( G4NucleiProperties::GetNuclearMass( static_cast<G4int>(theA+eps) , static_cast<G4int>(theZ+eps) )
	     ) / G4Neutron::Neutron()->GetPDGMass();
  
  G4ReactionProduct boosted;
  G4double aXsection;
  
  // MC integration loop
  G4int counter = 0;
  G4double buffer = 0;
  G4int size = G4int(std::max(10., aT/60*kelvin));
  G4ThreeVector neutronVelocity = 1./G4Neutron::Neutron()->GetPDGMass()*theNeutronRP.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();

  while(counter == 0 || std::abs(buffer-result/std::max(1,counter)) > 0.01*buffer) // Loop checking, 11.05.2015, T. Koi
  {
    if(counter) buffer = result/counter;
    while (counter<size) // Loop checking, 11.05.2015, T. Koi
    {
      counter ++;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus(eleMass, aT);
      boosted.Lorentz(theNeutronRP, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*((*theCrossSections)(index))).GetValue(theEkin, outOfRange);
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

G4int G4ParticleHPFissionData::GetVerboseLevel() const 
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}
void G4ParticleHPFissionData::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}
void G4ParticleHPFissionData::CrossSectionDescription(std::ostream& outFile) const
{
   outFile << "High Precision cross data based on Evaluated Nuclear Data Files (ENDF) for induced fission reaction of neutrons below 20MeV\n" ;
}
