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
// particle_hp -- source file
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
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleHPManager.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4ParticleHPData.hh"
#include "G4HadronicParameters.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"

G4ParticleHPInelasticData::G4ParticleHPInelasticData(G4ParticleDefinition* projectile)
  : G4VCrossSectionDataSet("")
{
  const char* dataDirVariable;
  G4String particleName;
  if( projectile == G4Neutron::Neutron() ) {
    dataDirVariable = "G4NEUTRONHPDATA";
  }else if( projectile == G4Proton::Proton() ) {
    dataDirVariable = "G4PROTONHPDATA";
    particleName = "Proton";
  }else if( projectile == G4Deuteron::Deuteron() ) {
    dataDirVariable = "G4DEUTERONHPDATA";
    particleName = "Deuteron";
  }else if( projectile == G4Triton::Triton() ) {
    dataDirVariable = "G4TRITONHPDATA";
    particleName = "Triton";
  }else if( projectile == G4He3::He3() ) {
    dataDirVariable = "G4HE3HPDATA";
    particleName = "He3";
  }else if( projectile == G4Alpha::Alpha() ) {
    dataDirVariable = "G4ALPHAHPDATA";
    particleName = "Alpha";
  } else {
    G4String message("G4ParticleHPInelasticData may only be called for neutron, proton, deuteron, triton, He3 or alpha, while it is called for " + projectile->GetParticleName());
    throw G4HadronicException(__FILE__, __LINE__,message.c_str());
  }
  G4String dataName = projectile->GetParticleName()+"HPInelasticXS";
  dataName.at(0) = (char)std::toupper(dataName.at(0)) ;
  SetName( dataName );

  if ( !G4FindDataDir(dataDirVariable) && !G4FindDataDir( "G4PARTICLEHPDATA" ) ){
    G4String message("Please setenv G4PARTICLEHPDATA (recommended) or, at least setenv " +
                     G4String(dataDirVariable) + " to point to the " + projectile->GetParticleName() + " cross-section files.");
    throw G4HadronicException(__FILE__, __LINE__,message.c_str());
  }

  G4String dirName;
  if ( G4FindDataDir(dataDirVariable) ) {
    dirName = G4FindDataDir(dataDirVariable);
  } else {
    G4String baseName = G4FindDataDir( "G4PARTICLEHPDATA" );
    dirName = baseName + "/" + particleName;
  }
  #ifdef G4VERBOSE
  if ( G4HadronicParameters::Instance()->GetVerboseLevel() > 0 ) {
    G4cout << "@@@ G4ParticleHPInelasticData instantiated for particle " << projectile->GetParticleName() << " data directory variable is " << dataDirVariable << " pointing to " << dirName << G4endl;
  }
  #endif

  SetMinKinEnergy( 0*CLHEP::MeV );                                   
  SetMaxKinEnergy( 20*CLHEP::MeV );                                   

  theCrossSections = 0;
  theProjectile=projectile;

  theHPData = nullptr;
  instanceOfWorker = false;
  if ( G4Threading::IsMasterThread() ) {
    theHPData = new G4ParticleHPData( theProjectile ); 
  } else {
    instanceOfWorker = true;
  }
  element_cache = nullptr;
  material_cache = nullptr;
  ke_cache = 0.0; 
  xs_cache = 0.0; 
}
   
G4ParticleHPInelasticData::~G4ParticleHPInelasticData()
{
   if ( theCrossSections != nullptr && instanceOfWorker != true ) {
     theCrossSections->clearAndDestroy();
     delete theCrossSections;
     theCrossSections = nullptr;
   }
   if ( theHPData != nullptr && instanceOfWorker != true ) {
     delete theHPData;
     theHPData = nullptr;
   }
}

G4bool G4ParticleHPInelasticData::IsIsoApplicable( const G4DynamicParticle* dp , 
                                                G4int /*Z*/ , G4int /*A*/ ,
                                                const G4Element* /*elm*/ ,
                                                const G4Material* /*mat*/ )
{
   G4double eKin = dp->GetKineticEnergy();
   if ( eKin > GetMaxKinEnergy() 
     || eKin < GetMinKinEnergy() 
     || dp->GetDefinition() != theProjectile ) return false;                                   

   return true;
}

G4double G4ParticleHPInelasticData::GetIsoCrossSection( const G4DynamicParticle* dp ,
                                   G4int /*Z*/ , G4int /*A*/ ,
                                   const G4Isotope* /*iso*/  ,
                                   const G4Element* element ,
                                   const G4Material* material )
{
   if ( dp->GetKineticEnergy() == ke_cache && element == element_cache
    &&  material == material_cache ) return xs_cache;

   ke_cache = dp->GetKineticEnergy();
   element_cache = element;
   material_cache = material;
   G4double xs = GetCrossSection( dp , element , material->GetTemperature() );
   xs_cache = xs;
   return xs;
}

void G4ParticleHPInelasticData::BuildPhysicsTable( const G4ParticleDefinition& projectile )
{
   if ( G4Threading::IsWorkerThread() ) {
      theCrossSections = G4ParticleHPManager::GetInstance()->GetInelasticCrossSections( &projectile );
      return;
   } else {
      if ( theHPData == nullptr ) theHPData = G4ParticleHPData::Instance( const_cast<G4ParticleDefinition*> ( &projectile ) ); 
   }

  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  if ( theCrossSections == nullptr ) 
    theCrossSections = new G4PhysicsTable( numberOfElements );
  else
    theCrossSections->clearAndDestroy();

  // make a PhysicsVector for each element

  static G4ThreadLocal G4ElementTable *theElementTable  = nullptr ;
  if (!theElementTable) theElementTable= G4Element::GetElementTable();
  for( std::size_t i=0; i<numberOfElements; ++i )
  {
    G4PhysicsVector* physVec = theHPData->MakePhysicsVector((*theElementTable)[i], this);
    theCrossSections->push_back(physVec);
  }
  G4ParticleHPManager::GetInstance()->RegisterInelasticCrossSections( &projectile , theCrossSections );
}

void G4ParticleHPInelasticData::DumpPhysicsTable(const G4ParticleDefinition& projectile)
{
  if(&projectile!=theProjectile) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use ParticleHP data for a wrong projectile!!!");  

#ifdef G4VERBOSE
  if ( G4HadronicParameters::Instance()->GetVerboseLevel() == 0 ) return;

  // Dump element based cross section
  // range 10e-5 eV to 20 MeV
  // 10 point per decade
  // in barn

  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "Inelastic Cross Section of Neutron HP"<< G4endl;
  G4cout << "(Pointwise cross-section at 0 Kelvin.)" << G4endl;
  G4cout << G4endl;
  G4cout << "Name of Element" << G4endl;
  G4cout << "Energy[eV]  XS[barn]" << G4endl;
  G4cout << G4endl;

  std::size_t numberOfElements = G4Element::GetNumberOfElements();
  static G4ThreadLocal G4ElementTable *theElementTable  = nullptr ;
  if (!theElementTable) theElementTable= G4Element::GetElementTable();

  for ( std::size_t i = 0 ; i < numberOfElements ; ++i )
  {
     G4cout << (*theElementTable)[i]->GetName() << G4endl;

     G4int ie = 0;

     for ( ie = 0 ; ie < 130 ; ie++ )
     {
       G4double eKinetic = 1.0e-5 * G4Pow::GetInstance()->powA ( 10.0 , ie/10.0 ) *CLHEP::eV;
       G4bool outOfRange = false;

       if ( eKinetic < 20*CLHEP::MeV )
       {
         G4cout << eKinetic/CLHEP::eV << " "
                << (*((*theCrossSections)(i))).GetValue(eKinetic, outOfRange)/CLHEP::barn
                << G4endl;
       }
     }
     G4cout << G4endl;
  }
#endif
}

G4double G4ParticleHPInelasticData::
GetCrossSection(const G4DynamicParticle* projectile, const G4Element*anE, G4double aT)
{
  G4double result = 0;
  G4bool outOfRange;
  std::size_t index = anE->GetIndex();

  // prepare neutron
  G4double eKinetic = projectile->GetKineticEnergy();

  if ( G4ParticleHPManager::GetInstance()->GetNeglectDoppler() )
  {
     //NEGLECT_DOPPLER
     G4double factor = 1.0;
     if ( eKinetic < aT * CLHEP::k_Boltzmann ) 
     {
        // below 0.1 eV neutrons 
        // Have to do some, but now just igonre.   
        // Will take care after performance check.  
        // factor = factor * targetV;
     }
     return ( (*((*theCrossSections)(index))).GetValue(eKinetic, outOfRange) )* factor;
  }   

  G4ReactionProduct theNeutron( projectile->GetDefinition() );
  theNeutron.SetMomentum( projectile->GetMomentum() );
  theNeutron.SetKineticEnergy( eKinetic );

  // prepare thermal nucleus
  G4Nucleus aNuc;
  G4double eps = 0.0001;
  G4double theA = anE->GetN();
  G4double theZ = anE->GetZ();
  G4double eleMass; 
  eleMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(theA+eps), static_cast<G4int>(theZ+eps) );
  
  G4ReactionProduct boosted;
  G4double aXsection;
  
  // MC integration loop
  G4int counter = 0;
  G4int failCount = 0;
  G4double buffer = 0; G4int size = G4int(std::max(10., aT/60*CLHEP::kelvin));
  G4ThreeVector neutronVelocity = 1./theProjectile->GetPDGMass()*theNeutron.GetMomentum();
  G4double neutronVMag = neutronVelocity.mag();

#ifndef G4PHP_DOPPLER_LOOP_ONCE
  while(counter == 0 || std::abs(buffer-result/std::max(1,counter)) > 0.01*buffer) // Loop checking, 11.05.2015, T. Koi
  {
    if(counter) buffer = result/counter;
    while (counter<size) // Loop checking, 11.05.2015, T. Koi
    {
      ++counter;
#endif
      G4ReactionProduct aThermalNuc = aNuc.GetThermalNucleus( eleMass/G4Neutron::Neutron()->GetPDGMass(), aT );
      boosted.Lorentz(theNeutron, aThermalNuc);
      G4double theEkin = boosted.GetKineticEnergy();
      aXsection = (*((*theCrossSections)(index))).GetValue(theEkin, outOfRange);
     if(aXsection <0) 
      {
        if(failCount<1000)
	{
	  ++failCount;
#ifndef G4PHP_DOPPLER_LOOP_ONCE
	  --counter;
	  continue;
#endif
	}
	else
	{
	  aXsection = 0;
	}
      }
      // velocity correction.
      G4ThreeVector targetVelocity = 1./aThermalNuc.GetMass()*aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity-neutronVelocity).mag()/neutronVMag;
      result += aXsection;
    }
#ifndef G4PHP_DOPPLER_LOOP_ONCE
    size += size;
  }
  result /= counter;
#endif

  return result;
}

G4int G4ParticleHPInelasticData::GetVerboseLevel() const 
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}

void G4ParticleHPInelasticData::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void G4ParticleHPInelasticData::CrossSectionDescription(std::ostream& outFile) const
{
   outFile << "Extension of High Precision cross section for inelastic reaction of proton, deuteron, triton, He3 and alpha below 20MeV\n";
}
