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
// Class Description
// Cross-section data set for a high precision (based on JENDL_HE evaluated data
// libraries) description of elastic scattering 20 MeV ~ 3 GeV;
// Class Description - End

// 15-Nov-06 First Implementation is done by T. Koi (SLAC/SCCS)
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPJENDLHEData.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ElementTable.hh"
#include "G4ParticleHPData.hh"
#include "G4Pow.hh"

G4bool G4ParticleHPJENDLHEData::IsApplicable(const G4DynamicParticle*aP, const G4Element* anE)
{
   G4bool result = true;
   G4double eKin = aP->GetKineticEnergy();
   //if(eKin>20*MeV||aP->GetDefinition()!=G4Neutron::Neutron()) result = false;
   if ( eKin < 20*MeV || 3*GeV < eKin || aP->GetDefinition()!=G4Neutron::Neutron() ) 
   {
      result = false;
   } 
// Element Check 
   else if ( !(vElement[ anE->GetIndex() ]) ) result = false;

   return result;
}


G4ParticleHPJENDLHEData::G4ParticleHPJENDLHEData()
{
  for ( std::map< G4int , std::map< G4int , G4PhysicsVector* >* >::iterator itZ = mIsotope.begin();
        itZ != mIsotope.end(); ++itZ ) {
    std::map< G4int , G4PhysicsVector* >* pointer_map = itZ->second;
    if ( pointer_map ) {
      for ( std::map< G4int , G4PhysicsVector* >::iterator itA = pointer_map->begin();
            itA != pointer_map->end() ; ++itA ) {
        G4PhysicsVector* pointerPhysicsVector = itA->second;
        if ( pointerPhysicsVector ) {
          delete pointerPhysicsVector;
          itA->second = NULL;
        }
      }
      delete pointer_map;
      itZ->second = NULL;
    }
  }
  mIsotope.clear();
}


G4ParticleHPJENDLHEData::G4ParticleHPJENDLHEData( G4String reaction , G4ParticleDefinition* pd )
  : G4VCrossSectionDataSet( "JENDLHE"+reaction+"CrossSection" )
{
   reactionName = reaction;
   BuildPhysicsTable( *pd );
}


G4ParticleHPJENDLHEData::~G4ParticleHPJENDLHEData()
{
}


void G4ParticleHPJENDLHEData::BuildPhysicsTable( const G4ParticleDefinition& aP )
{
   particleName = aP.GetParticleName();

   G4String baseName = G4FindDataDir( "G4NEUTRONHPDATA" );
   G4String dirName = baseName+"/JENDL_HE/"+particleName+"/"+reactionName ;
   G4String aFSType = "/CrossSection/";
   G4ParticleHPNames theNames; 

   G4String filename;

   // Create JENDL_HE data 
   // Create map element or isotope  

   std::size_t numberOfElements = G4Element::GetNumberOfElements();

   // make a PhysicsVector for each element

   static G4ThreadLocal G4ElementTable *theElementTable  = 0 ; if (!theElementTable) theElementTable= G4Element::GetElementTable();
   vElement.clear();
   vElement.resize( numberOfElements );
   for ( std::size_t i = 0; i < numberOfElements; ++i )
   {
      G4Element* theElement = (*theElementTable)[i];
      vElement[i] = false;

      // isotope
      G4int nIso = (G4int)(*theElementTable)[i]->GetNumberOfIsotopes();
      G4int Z = (G4int)(*theElementTable)[i]->GetZ();
      if ( nIso!=0 )
      {
         G4bool found_at_least_one = false; 
         for ( G4int i1 = 0; i1 < nIso; ++i1 )
         {
             G4int A = theElement->GetIsotope(i1)->GetN();

             if ( isThisNewIsotope( Z , A ) ) 
             {
                std::stringstream ss; 
                ss << dirName << aFSType << Z << "_" << A << "_" << theNames.GetName( Z-1 );
                filename = ss.str();
                std::fstream file;
                file.open ( filename , std::fstream::in );
                G4int dummy;
                file >> dummy;
                if ( file.good() ) 
                {
                   found_at_least_one = true;

                   // read the file
                   G4PhysicsVector* aPhysVec = readAFile ( &file );
                   registAPhysicsVector( Z , A , aPhysVec );
                }
                else 
                { 
                   //G4cout << "No file for "<< reactionType << " Z=" << Z << ", A=" << A << G4endl;
                }
                file.close();
             }
             else
             {
                found_at_least_one = TRUE;
             }
          }
          if ( found_at_least_one ) vElement[i] = true;
       }
       else
       {
          G4StableIsotopes theStableOnes;
          G4int first = theStableOnes.GetFirstIsotope( Z );
          G4bool found_at_least_one = FALSE; 
          for ( G4int i1 = 0; i1 < theStableOnes.GetNumberOfIsotopes( static_cast<G4int>(theElement->GetZ() ) ); i1++)
          {
             G4int A = theStableOnes.GetIsotopeNucleonCount( first+i1 );
             if ( isThisNewIsotope( Z , A ) ) 
             {

                std::stringstream ss; 
                ss << dirName << aFSType << Z << "_" << A << "_" << theNames.GetName( Z-1 );
                filename = ss.str();

                std::fstream file;
                file.open ( filename , std::fstream::in );
                G4int dummy;
                file >> dummy;
                if ( file.good() ) 
                {
                   //G4cout << "Found file for Z=" << Z << ", A=" << A << ", as " << filename << G4endl;
                   found_at_least_one = TRUE;
                   //Read the file

                   G4PhysicsVector* aPhysVec = readAFile ( &file );

                   //Regist the PhysicsVector
                   registAPhysicsVector( Z , A , aPhysVec );
                }
                else 
                { 
                   //G4cout << "No file for "<< reactionType << " Z=" << Z << ", A=" << A << G4endl;
                }
                file.close();
             }
             else
             {
                found_at_least_one = TRUE;
             }
          }

          if ( found_at_least_one ) vElement[i] = true;
       }
   }
}


void G4ParticleHPJENDLHEData::DumpPhysicsTable(const G4ParticleDefinition& aP)
{
  if(&aP!=G4Neutron::Neutron()) 
     throw G4HadronicException(__FILE__, __LINE__, "Attempt to use NeutronHP data for particles other than neutrons!!!");  
}


G4double G4ParticleHPJENDLHEData::
GetCrossSection(const G4DynamicParticle* aP, const G4Element*anE, G4double )
{
   // Primary energy >20MeV
   // Thus not taking into account of Doppler broadening 
   // also not taking into account of Target thermal motions

   G4double result = 0;

   G4double ek = aP->GetKineticEnergy();

   G4int nIso = (G4int)anE->GetNumberOfIsotopes();
   G4int Z = (G4int)anE->GetZ();
   if ( nIso!=0 )
   {
      for ( G4int i1 = 0; i1 < nIso; ++i1 )
      {
         G4int A = anE->GetIsotope(i1)->GetN();
         G4double frac = anE->GetRelativeAbundanceVector()[ i1 ]; // This case does NOT request "*perCent".
         result += frac * getXSfromThisIsotope( Z , A , ek );
      }
   }
   else
   {
      G4StableIsotopes theStableOnes;
      G4int first = theStableOnes.GetFirstIsotope( Z );
      for ( G4int i1 = 0; i1 < theStableOnes.GetNumberOfIsotopes( (G4int)anE->GetZ() ); ++i1)
      {
         G4int A = theStableOnes.GetIsotopeNucleonCount( first+i1 );
         G4double frac = theStableOnes.GetAbundance( first+i1 )*perCent; // This case requests "*perCent".
         result += frac * getXSfromThisIsotope( Z , A , ek );
      }
   }
   return result;
}


G4PhysicsVector* G4ParticleHPJENDLHEData::readAFile ( std::fstream* file )
{
   G4int dummy;
   G4int len;
   *file >> dummy;
   *file >> len;

   std::vector< G4double > v_e; 
   std::vector< G4double > v_xs; 

   for ( G4int i = 0 ; i < len ; ++i )
   {
      G4double e;
      G4double xs;

      *file >> e; 
      *file >> xs;
      // data are written in eV and barn.    
      v_e.push_back( e*eV );
      v_xs.push_back( xs*barn );
   }

   G4PhysicsFreeVector* aPhysVec = new G4PhysicsFreeVector( static_cast< std::size_t >( len ) , v_e.front() , v_e.back() );

   for ( G4int i = 0 ; i < len ; ++i )
   {
      aPhysVec->PutValues( static_cast< std::size_t >( i ) , v_e[ i ] , v_xs[ i ] );   
   }

   return aPhysVec;
}


G4bool G4ParticleHPJENDLHEData::isThisInMap( G4int z , G4int a )
{
   if ( mIsotope.find ( z ) == mIsotope.end() ) return false;
   if ( mIsotope.find ( z ) -> second->find ( a ) ==  mIsotope.find ( z ) -> second->end() ) return false;
   return true; 
}


void G4ParticleHPJENDLHEData::registAPhysicsVector( G4int Z , G4int A , G4PhysicsVector* aPhysVec )
{
    std::pair< G4int , G4PhysicsVector* > aPair = std::pair < G4int , G4PhysicsVector* > ( A , aPhysVec );
    auto itm = mIsotope.find ( Z );
    if ( itm !=  mIsotope.cend() ) 
    { 
       itm->second->insert ( aPair ); 
    }  
    else
    {
       std::map< G4int , G4PhysicsVector* >* aMap = new std::map< G4int , G4PhysicsVector* >;
       aMap->insert ( aPair ); 
       mIsotope.insert( std::pair< G4int , std::map< G4int , G4PhysicsVector* >* > ( Z , aMap ) );
    }
}


G4double G4ParticleHPJENDLHEData::getXSfromThisIsotope( G4int Z , G4int A , G4double ek )
{
   G4double aXSection = 0.0;
   G4bool outOfRange;

   G4PhysicsVector* aPhysVec;
   if ( mIsotope.find ( Z )->second->find ( A ) != mIsotope.find ( Z )->second->end() )
   {
      aPhysVec = mIsotope.find ( Z )->second->find ( A )->second;
      aXSection = aPhysVec->GetValue( ek , outOfRange );
   } 
   else
   {
      // Select closest one in the same Z
      G4int delta0 = 99; // no mean for 99 
      for ( auto it = mIsotope.find ( Z )->second->cbegin();
                 it!= mIsotope.find ( Z )->second->cend(); ++it )
      {
         G4int delta = std::abs( A - it->first );
         if ( delta < delta0 ) delta0 = delta;
      }

      // Randomize of selection larger or smaller than A
      if ( G4UniformRand() < 0.5 ) delta0 *= -1;
      G4int A1 = A + delta0;
      if ( mIsotope.find ( Z )->second->find ( A1 ) != mIsotope.find ( Z )->second->cend() )
      {
         aPhysVec = mIsotope.find ( Z )->second->find ( A1 )->second;
      }
      else
      {
         A1 = A - delta0;
         aPhysVec = mIsotope.find ( Z )->second->find ( A1 )->second;
      }

      aXSection = aPhysVec->GetValue( ek , outOfRange );
      // X^(2/3) factor
      aXSection *= G4Pow::GetInstance()->A23( 1.0*A/ A1 );
   }

   return aXSection;
}
