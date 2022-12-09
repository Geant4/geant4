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
// 070523 add neglecting doppler broadening on the fly. T. Koi
// 070613 fix memory leaking by T. Koi
// 071002 enable cross section dump by T. Koi
// 080428 change checking point of "neglecting doppler broadening" flag 
//        from GetCrossSection to BuildPhysicsTable by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4ParticleHPData.hh"
#include "G4ParticleHPManager.hh"
#include "G4Threading.hh"
#include "G4HadronicParameters.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"

G4ParticleHPCaptureData::G4ParticleHPCaptureData()
:G4VCrossSectionDataSet("NeutronHPCaptureXS")
{
   SetMinKinEnergy( 0*MeV );                                   
   SetMaxKinEnergy( 20*MeV );                                   

   theCrossSections = 0;

   instanceOfWorker = false;
   if ( G4Threading::IsWorkerThread() ) {
      instanceOfWorker = true;
   }

   element_cache = nullptr;
   material_cache = nullptr;
   ke_cache = 0.0; 
   xs_cache = 0.0; 
    
   //BuildPhysicsTable(*G4Neutron::Neutron());
}
   
G4ParticleHPCaptureData::~G4ParticleHPCaptureData()
{
   if ( theCrossSections != nullptr && instanceOfWorker != true ) {
     theCrossSections->clearAndDestroy();
     delete theCrossSections;
     theCrossSections = nullptr;
   }
}
   
G4bool G4ParticleHPCaptureData::IsIsoApplicable( const G4DynamicParticle* dp , 
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

G4double G4ParticleHPCaptureData::GetIsoCrossSection( const G4DynamicParticle* dp ,
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
G4bool G4ParticleHPCaptureData::IsApplicable(const G4DynamicParticle*aP, const G4Element*)
{
  G4bool result = true;
  G4double eKin = aP->GetKineticEnergy();
  if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
  return result;
}
*/

void G4ParticleHPCaptureData::BuildPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

   if ( G4Threading::IsWorkerThread() ) {
      theCrossSections = G4ParticleHPManager::GetInstance()->GetCaptureCrossSections();
      return;
   }
  
  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  // G4cout << "CALLED G4ParticleHPCaptureData::BuildPhysicsTable "<<numberOfElements<<G4endl;
   // TKDB
   //if ( theCrossSections == 0 ) theCrossSections = new G4PhysicsTable( numberOfElements );
   if ( theCrossSections == nullptr ) 
      theCrossSections = new G4PhysicsTable( numberOfElements );
   else
      theCrossSections->clearAndDestroy();

  // make a PhysicsVector for each element

  static G4ThreadLocal G4ElementTable *theElementTable  = 0 ; if (!theElementTable) theElementTable= G4Element::GetElementTable();
  for( std::size_t i=0; i<numberOfElements; ++i )
  {
     #ifdef G4VERBOSE
     if(std::getenv("CaptureDataIndexDebug"))
     {
       std::size_t index_debug = ((*theElementTable)[i])->GetIndex();
       if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 0 ) G4cout << "IndexDebug "<< i <<" "<<index_debug<<G4endl;
     }
     #endif
     G4PhysicsVector* physVec = G4ParticleHPData::
      Instance(G4Neutron::Neutron())->MakePhysicsVector((*theElementTable)[i], this);
     theCrossSections->push_back(physVec);
  }

  G4ParticleHPManager::GetInstance()->RegisterCaptureCrossSections( theCrossSections );
}

void G4ParticleHPCaptureData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");

  #ifdef G4VERBOSE
  if ( G4HadronicParameters::Instance()->GetVerboseLevel() == 0 ) return;
  
//
// Dump element based cross section
// range 10e-5 eV to 20 MeV
// 10 point per decade
// in barn
//

   G4cout << G4endl;
   G4cout << G4endl;
   G4cout << "Capture Cross Section of Neutron HP"<< G4endl;
   G4cout << "(Pointwise cross-section at 0 Kelvin.)" << G4endl;
   G4cout << G4endl;
   G4cout << "Name of Element" << G4endl;
   G4cout << "Energy[eV]  XS[barn]" << G4endl;
   G4cout << G4endl;

   std::size_t numberOfElements = G4Element::GetNumberOfElements();
   static G4ThreadLocal G4ElementTable *theElementTable  = 0 ;
   if (!theElementTable) theElementTable= G4Element::GetElementTable();

   for ( std::size_t i = 0 ; i < numberOfElements ; ++i )
   {

      G4cout << (*theElementTable)[i]->GetName() << G4endl;

      G4int ie = 0;

      for ( ie = 0 ; ie < 130 ; ++ie )
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

   //G4cout << "G4ParticleHPCaptureData::DumpPhysicsTable still to be implemented"<<G4endl;
   #endif
}

G4double G4ParticleHPCaptureData::
GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double aT)
{
  G4double result = 0;
  G4bool outOfRange;
  G4int index = (G4int)anE->GetIndex();

  // prepare neutron
  G4double eKinetic = aP->GetKineticEnergy();

  if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() )
  {
     //NEGLECT_DOPPLER
     G4double factor = 1.0;
     if ( eKinetic < aT * k_Boltzmann ) 
     {
        // below 0.1 eV neutrons 
        // Have to do some, but now just igonre.   
        // Will take care after performance check.  
        // factor = factor * targetV;
     }
     return ( (*((*theCrossSections)(index))).GetValue(eKinetic, outOfRange) )* factor; 
  }

  G4ReactionProduct theNeutron( aP->GetDefinition() );
  theNeutron.SetMomentum( aP->GetMomentum() );
  theNeutron.SetKineticEnergy( eKinetic );

  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eps = 0.0001;
  G4double theA = anE->GetN();
  G4double theZ = anE->GetZ();
  G4double eleMass; 
  eleMass = G4NucleiProperties::GetNuclearMass( static_cast<G4int>(theA+eps) , static_cast<G4int>(theZ+eps) ) / G4Neutron::Neutron()->GetPDGMass();
  
  G4ReactionProduct boosted;
  G4double aXsection;
  
  // MC integration loop
  G4int counter = 0;
  G4double buffer = 0;
  G4int size = G4int(std::max(10., aT/60*kelvin));
  G4ThreeVector neutronVelocity = 1./G4Neutron::Neutron()->GetPDGMass()*theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();

  while(counter == 0 || std::abs(buffer-result/std::max(1,counter)) > 0.03*buffer) // Loop checking, 11.05.2015, T. Koi
  {
    if(counter) buffer = result/counter;
    while (counter<size) // Loop checking, 11.05.2015, T. Koi
    {
      ++counter;
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus(eleMass, aT);
      boosted.Lorentz(theNeutron, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*((*theCrossSections)(index))).GetValue(theEkin, outOfRange);
      // velocity correction, or luminosity factor...
      G4ThreeVector targetVelocity = 1./aThermalNuc.GetMass()*aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity-neutronVelocity).mag()/neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
/*
  // Checking impact of  G4NEUTRONHP_NEGLECT_DOPPLER
  G4cout << " result " << result << " " 
         << (*((*theCrossSections)(index))).GetValue(eKinetic, outOfRange) << " " 
         << (*((*theCrossSections)(index))).GetValue(eKinetic, outOfRange) /result << G4endl;
*/
  return result;
}

G4int G4ParticleHPCaptureData::GetVerboseLevel() const 
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}

void G4ParticleHPCaptureData::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void G4ParticleHPCaptureData::CrossSectionDescription(std::ostream& outFile) const
{
   outFile << "High Precision cross data based on Evaluated Nuclear Data Files (ENDF) for radiative capture reaction of neutrons below 20MeV\n" ;
}
