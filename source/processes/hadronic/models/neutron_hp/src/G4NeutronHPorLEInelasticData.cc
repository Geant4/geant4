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

#include "G4NeutronHPorLEInelasticData.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4NeutronHPData.hh"

#include "G4PhysicsVector.hh"



G4NeutronHPorLEInelasticData::G4NeutronHPorLEInelasticData( G4NeutronHPChannelList* pChannel , std::set< G4String >* pSet )
{
   theInelasticChannel = pChannel;
   unavailable_elements = pSet;   
   BuildPhysicsTable(*G4Neutron::Neutron());
}

 

G4bool G4NeutronHPorLEInelasticData::IsApplicable(const G4DynamicParticle*aP, const G4Element* anElement)
{
   G4bool result = true;
   G4double eKin = aP->GetKineticEnergy();
   if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
   if ( unavailable_elements->find( anElement->GetName() ) != unavailable_elements->end() ) result = false;
   return result;
}

G4NeutronHPorLEInelasticData::G4NeutronHPorLEInelasticData()
{
//   BuildPhysicsTable(*G4Neutron::Neutron());
}

 

G4NeutronHPorLEInelasticData::~G4NeutronHPorLEInelasticData()
{
//  delete theCrossSections;
}


#include "G4NeutronHPInelasticData.hh"
#include "G4LPhysicsFreeVector.hh"
//#include "G4NeutronHPElementData.hh"
   
void G4NeutronHPorLEInelasticData::BuildPhysicsTable( const G4ParticleDefinition& aP )
{
   if( &aP!=G4Neutron::Neutron() ) 
      throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  

   size_t numberOfElements = G4Element::GetNumberOfElements();
   theCrossSections = new G4PhysicsTable( numberOfElements );

   static const G4ElementTable *theElementTable = G4Element::GetElementTable(); 
   for ( size_t i=0 ; i < numberOfElements; ++i )
   {
      G4PhysicsVector* thePhysVec = new G4LPhysicsFreeVector(0, 0, 0);

      if ( unavailable_elements->find( (*theElementTable)[i]->GetName() ) == unavailable_elements->end() ) 
      { 

         G4NeutronHPElementData* theElementData = new G4NeutronHPElementData();
         theElementData->Init( (*theElementTable)[i] );  

         G4NeutronHPVector* theHPVector = theElementData->GetData( (G4NeutronHPInelasticData*)this ); 

         G4int len = theHPVector->GetVectorLength();

         if ( len!=0 ) 
         {
            G4double emin = theHPVector->GetX(0);
            G4double emax = theHPVector->GetX(len-1);

            G4LPhysicsFreeVector* aPhysVector= new G4LPhysicsFreeVector ( len , emin , emax );
            for ( G4int i=0; i<len; i++ )
            {
               aPhysVector->PutValues( i , theHPVector->GetX(i) , theHPVector->GetY(i) );
            }
            delete thePhysVec;
            thePhysVec = aPhysVector;
         }

         //G4PhysicsVector* physVec = G4NeutronHPData::
         //Instance()->MakePhysicsVector((*theElementTable)[i], this);
         //theCrossSections->push_back(physVec);
      }

      theCrossSections->push_back(thePhysVec);
   }
}
 


void G4NeutronHPorLEInelasticData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
//  G4cout << "G4NeutronHPorLEInelasticData::DumpPhysicsTable still to be implemented"<<G4endl;
}



#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"

G4double G4NeutronHPorLEInelasticData::
GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double aT)
{

  // G4cout << "Choice G4NeutronHPorLEInelasticData for element " << anE->GetName() << G4endl;
  G4double result = 0;
  //G4bool outOfRange;
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
      aXsection = theInelasticChannel[index].GetXsec( theEkin );
      // velocity correction.
      G4ThreeVector targetVelocity = 1./aThermalNuc.GetMass()*aThermalNuc.GetMomentum();
      aXsection *= (targetVelocity-neutronVelocity).mag()/neutronVMag;
      result += aXsection;
    }
    size += size;
  }
  result /= counter;
  //return result;
  return result*barn;
}
