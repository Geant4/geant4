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
// Thermal Neutron Scattering
// Koi, Tatsumi (SLAC/SCCS)
//
// Class Description:
//
// Final State Generators for a high precision (based on evaluated data
// libraries) description of themal neutron scattering below 4 eV;
// Based on Thermal neutron scattering files
// from the evaluated nuclear data files ENDF/B-VI, Release2
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.


// 070625 Fix memory leaking at destructor by T. Koi 
// 081201 Fix memory leaking at destructor by T. Koi 
// 100729 Add model name in constructor Problem #1116
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPThermalScattering.hh"
#include "G4ParticleHPThermalScatteringData.hh"
#include "G4ParticleHPThermalScatteringNames.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"
#include "G4Threading.hh"

G4ParticleHPThermalScattering::G4ParticleHPThermalScattering()
                             :G4HadronicInteraction("NeutronHPThermalScattering")
,coherentFSs(NULL)
,incoherentFSs(NULL)
,inelasticFSs(NULL)
{
   theHPElastic = new G4ParticleHPElastic();

   SetMinEnergy( 0.*eV );
   SetMaxEnergy( 4*eV );
   theXSection = new G4ParticleHPThermalScatteringData();

   nMaterial = 0;
   nElement = 0;
}



G4ParticleHPThermalScattering::~G4ParticleHPThermalScattering()
{
   delete theHPElastic;
}



void G4ParticleHPThermalScattering::clearCurrentFSData() {

if ( incoherentFSs != NULL ) {
   for ( std::map < G4int , std::map < G4double , std::vector < E_isoAng* >* >* >::iterator it = incoherentFSs->begin() ; it != incoherentFSs->end() ; it++ )
   {
      std::map < G4double , std::vector < E_isoAng* >* >::iterator itt;
      for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
      {
         std::vector< E_isoAng* >::iterator ittt;
         for ( ittt = itt->second->begin(); ittt != itt->second->end() ; ittt++ )
         {
            delete *ittt;
         }
         delete itt->second;
      }
      delete it->second;
   }
}

if ( coherentFSs != NULL ) {
   for ( std::map < G4int , std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* >::iterator it = coherentFSs->begin() ; it != coherentFSs->end() ; it++ )
   {
      std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >::iterator itt;
      for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
      {
         std::vector < std::pair< G4double , G4double >* >::iterator ittt;
         for ( ittt = itt->second->begin(); ittt != itt->second->end() ; ittt++ )
         {
            delete *ittt;
         }
         delete itt->second;
      }
      delete it->second;
   }
}

if ( inelasticFSs != NULL ) {
   for ( std::map < G4int ,  std::map < G4double , std::vector < E_P_E_isoAng* >* >* >::iterator it = inelasticFSs->begin() ; it != inelasticFSs->end() ; it++ )
   {
      std::map < G4double , std::vector < E_P_E_isoAng* >* >::iterator itt;
      for ( itt = it->second->begin() ; itt != it->second->end() ; itt++ )
      {
         std::vector < E_P_E_isoAng* >::iterator ittt;
         for ( ittt = itt->second->begin(); ittt != itt->second->end() ; ittt++ )
         {
            std::vector < E_isoAng* >::iterator it4;
            for ( it4 = (*ittt)->vE_isoAngle.begin() ; it4 != (*ittt)->vE_isoAngle.end() ; it4++ )
            {
               delete *it4;
            }
            delete *ittt;
         }
         delete itt->second;
      }
      delete it->second;
   }
}

   incoherentFSs = NULL;
   coherentFSs = NULL;
   inelasticFSs = NULL;

}



void G4ParticleHPThermalScattering::BuildPhysicsTable(const G4ParticleDefinition& particle) {
   buildPhysicsTable();
   theHPElastic->BuildPhysicsTable( particle );
}



std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* G4ParticleHPThermalScattering::readACoherentFSDATA( G4String name )
{

   std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* aCoherentFSDATA = new std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >;

   //std::ifstream theChannel( name.c_str() );
   std::istringstream theChannel(std::ios::in);
   G4ParticleHPManager::GetInstance()->GetDataStream(name,theChannel);

   std::vector< G4double > vBraggE;

   G4int dummy; 
   while ( theChannel >> dummy )   // MF // Loop checking, 11.05.2015, T. Koi
   {
      theChannel >> dummy;   // MT
      G4double temp; 
      theChannel >> temp;   
      std::vector < std::pair< G4double , G4double >* >* anBragE_P = new std::vector < std::pair< G4double , G4double >* >;
     
      G4int n; 
      theChannel >> n;   
      for ( G4int i = 0 ; i < n ; i++ )
      {
          G4double Ei; 
          G4double Pi;
          if ( aCoherentFSDATA->size() == 0 ) 
          {
             theChannel >> Ei;
             vBraggE.push_back( Ei );
          } 
          else 
          {
             Ei = vBraggE[ i ]; 
          } 
          theChannel >> Pi;   
          anBragE_P->push_back ( new std::pair < G4double , G4double > ( Ei , Pi ) );
          //G4cout << "Coherent Elastic " << Ei << " " << Pi << G4endl;   
      }
      aCoherentFSDATA->insert ( std::pair < G4double , std::vector < std::pair< G4double , G4double >* >*  > ( temp , anBragE_P ) );
   }
 
   return aCoherentFSDATA;
}



std::map < G4double , std::vector < E_P_E_isoAng* >* >* G4ParticleHPThermalScattering::readAnInelasticFSDATA ( G4String name )
{
   std::map < G4double , std::vector < E_P_E_isoAng* >* >* anT_E_P_E_isoAng = new std::map < G4double , std::vector < E_P_E_isoAng* >* >;

   //std::ifstream theChannel( name.c_str() );
   std::istringstream theChannel(std::ios::in);
   G4ParticleHPManager::GetInstance()->GetDataStream(name,theChannel);

   G4int dummy; 
   while ( theChannel >> dummy )   // MF // Loop checking, 11.05.2015, T. Koi
   {
      theChannel >> dummy;   // MT
      G4double temp; 
      theChannel >> temp;   
      std::vector < E_P_E_isoAng* >* vE_P_E_isoAng = new std::vector < E_P_E_isoAng* >;
      G4int n;
      theChannel >> n;   
      for ( G4int i = 0 ; i < n ; i++ )
      {
          vE_P_E_isoAng->push_back ( readAnE_P_E_isoAng ( &theChannel ) );
      }
      anT_E_P_E_isoAng->insert ( std::pair < G4double , std::vector < E_P_E_isoAng* >* > ( temp , vE_P_E_isoAng ) );
   }    
   //theChannel.close();

   return anT_E_P_E_isoAng; 
}



E_P_E_isoAng* G4ParticleHPThermalScattering::readAnE_P_E_isoAng( std::istream* file )
{
   E_P_E_isoAng* aData = new E_P_E_isoAng;

   G4double dummy;
   G4double energy;
   G4int nep , nl;
   *file >> dummy;
   *file >> energy;
   aData->energy = energy*eV;
   *file >> dummy;
   *file >> dummy;
   *file >> nep;
   *file >> nl;
   aData->n = nep/nl;
   for ( G4int i = 0 ; i < aData->n ; i++ )
   {
      G4double prob;
      E_isoAng* anE_isoAng = new E_isoAng;
      aData->vE_isoAngle.push_back( anE_isoAng );
      *file >> energy;
      anE_isoAng->energy = energy*eV; 
      anE_isoAng->n = nl - 2;  
      anE_isoAng->isoAngle.resize( anE_isoAng->n ); 
      *file >> prob;
      aData->prob.push_back( prob );
      //G4cout << "G4ParticleHPThermalScattering inelastic " << energy/eV << " " <<  i << " " << prob << " " << aData->prob[ i ] << G4endl; 
      for ( G4int j = 0 ; j < anE_isoAng->n ; j++ )
      {
         G4double x;
         *file >> x;
         anE_isoAng->isoAngle[j] = x ;
         //G4cout << "G4ParticleHPThermalScattering inelastic " << x << anE_isoAng->isoAngle[j] << G4endl; 
      }
   } 

   // Calcuate sum_of_provXdEs
   G4double total = 0;  
   for ( G4int i = 0 ; i < aData->n - 1 ; i++ )
   {
      G4double E_L = aData->vE_isoAngle[i]->energy/eV;
      G4double E_H = aData->vE_isoAngle[i+1]->energy/eV;
      G4double dE = E_H - E_L;
      total += ( ( aData->prob[i] ) * dE );
   }
   aData->sum_of_probXdEs = total;

   return aData;
}



std::map < G4double , std::vector < E_isoAng* >* >* G4ParticleHPThermalScattering::readAnIncoherentFSDATA ( G4String name )
{
   std::map < G4double , std::vector < E_isoAng* >* >* T_E = new std::map < G4double , std::vector < E_isoAng* >* >;

   //std::ifstream theChannel( name.c_str() );
   std::istringstream theChannel(std::ios::in);
   G4ParticleHPManager::GetInstance()->GetDataStream(name,theChannel);

   G4int dummy; 
   while ( theChannel >> dummy )   // MF // Loop checking, 11.05.2015, T. Koi
   {
      theChannel >> dummy;   // MT
      G4double temp; 
      theChannel >> temp;   
      std::vector < E_isoAng* >* vE_isoAng = new std::vector < E_isoAng* >;
      G4int n;
      theChannel >> n;   
      for ( G4int i = 0 ; i < n ; i++ )
        vE_isoAng->push_back ( readAnE_isoAng( &theChannel ) );
      T_E->insert ( std::pair < G4double , std::vector < E_isoAng* >* > ( temp , vE_isoAng ) );
   }
   //theChannel.close();

   return T_E;
}



E_isoAng* G4ParticleHPThermalScattering::readAnE_isoAng( std::istream* file )
{
   E_isoAng* aData = new E_isoAng;

   G4double dummy;
   G4double energy;
   G4int n;
   *file >> dummy;
   *file >> energy;
   *file >> dummy;
   *file >> dummy;
   *file >> n;
   *file >> dummy;
   aData->energy = energy*eV;
   aData->n = n-2;
   aData->isoAngle.resize( n );

   *file >> dummy;
   *file >> dummy;
   for ( G4int i = 0 ; i < aData->n ; i++ )
      *file >> aData->isoAngle[i];

   return aData;
}



G4HadFinalState* G4ParticleHPThermalScattering::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
{

// Select Element > Reaction >

   const G4Material * theMaterial = aTrack.GetMaterial();
   G4double aTemp = theMaterial->GetTemperature();
   G4int n = theMaterial->GetNumberOfElements();

   G4bool findThermalElement = false;
   G4int ielement;
   const G4Element* theElement = NULL;
   for ( G4int i = 0; i < n ; i++ )
   {
      theElement = theMaterial->GetElement(i);
      //Select target element 
      if ( aNucleus.GetZ_asInt() == (G4int)(theElement->GetZ() + 0.5 ) )
      {
         //Check Applicability of Thermal Scattering 
         if (  getTS_ID( NULL , theElement ) != -1 )
         {
            ielement = getTS_ID( NULL , theElement );
            findThermalElement = true;
            break;
         }
         else if (  getTS_ID( theMaterial , theElement ) != -1 )
         {
            ielement = getTS_ID( theMaterial , theElement );
            findThermalElement = true;
            break;
         }
      }       
   } 

   if ( findThermalElement == true )
   {

//    Select Reaction  (Inelastic, coherent, incoherent)  

      const G4ParticleDefinition* pd = aTrack.GetDefinition();
      G4DynamicParticle* dp = new G4DynamicParticle ( pd , aTrack.Get4Momentum() );
      G4double total = theXSection->GetCrossSection( dp , theElement , theMaterial );
      G4double inelastic = theXSection->GetInelasticCrossSection( dp , theElement , theMaterial );
   

      G4double random = G4UniformRand();
      if ( random <= inelastic/total ) 
      {
         // Inelastic

         // T_L and T_H 
         std::map < G4double , std::vector< E_P_E_isoAng* >* >::iterator it; 
         std::vector<G4double> v_temp;
         v_temp.clear();
         for ( it = inelasticFSs->find( ielement )->second->begin() ; it != inelasticFSs->find( ielement )->second->end() ; it++ )
         {
            v_temp.push_back( it->first );
         }

//                   T_L         T_H 
         std::pair < G4double , G4double > tempLH = find_LH ( aTemp , &v_temp );
//
//       For T_L aNEP_EPM_TL  and T_H aNEP_EPM_TH
//
         std::vector< E_P_E_isoAng* >* vNEP_EPM_TL = 0;
         std::vector< E_P_E_isoAng* >* vNEP_EPM_TH = 0;

         if ( tempLH.first != 0.0 && tempLH.second != 0.0 ) 
         {
            vNEP_EPM_TL = inelasticFSs->find( ielement )->second->find ( tempLH.first/kelvin )->second;
            vNEP_EPM_TH = inelasticFSs->find( ielement )->second->find ( tempLH.second/kelvin )->second;
         }
         else if ( tempLH.first == 0.0 )
         {
            std::map < G4double , std::vector< E_P_E_isoAng* >* >::iterator itm;  
            itm = inelasticFSs->find( ielement )->second->begin();
            vNEP_EPM_TL = itm->second;
            itm++;
            vNEP_EPM_TH = itm->second;
            tempLH.first = tempLH.second;
            tempLH.second = itm->first;
         }
         else if (  tempLH.second == 0.0 )
         {
            std::map < G4double , std::vector< E_P_E_isoAng* >* >::iterator itm;  
            itm = inelasticFSs->find( ielement )->second->end();
            itm--;
            vNEP_EPM_TH = itm->second;
            itm--;
            vNEP_EPM_TL = itm->second;
            tempLH.second = tempLH.first;
            tempLH.first = itm->first;
         } 

         G4double rand_for_sE = G4UniformRand();

         std::pair< G4double , E_isoAng > TL = create_sE_and_EPM_from_pE_and_vE_P_E_isoAng ( rand_for_sE ,  aTrack.GetKineticEnergy() , vNEP_EPM_TL );
         std::pair< G4double , E_isoAng > TH = create_sE_and_EPM_from_pE_and_vE_P_E_isoAng ( rand_for_sE ,  aTrack.GetKineticEnergy() , vNEP_EPM_TH );

         G4double sE;
         sE = get_linear_interpolated ( aTemp , std::pair < G4double , G4double > ( tempLH.first , TL.first ) , std::pair < G4double , G4double > ( tempLH.second , TH.first ) );  

         G4double mu=1.0;
         E_isoAng anE_isoAng; 
         if ( TL.second.n == TH.second.n ) 
         {
            anE_isoAng.energy = sE; 
            anE_isoAng.n =  TL.second.n; 
            for ( G4int i=0 ; i < anE_isoAng.n ; i++ )
            { 
               G4double angle;
               angle = get_linear_interpolated ( aTemp , std::pair< G4double , G4double > (  tempLH.first , TL.second.isoAngle[ i ] ) , std::pair< G4double , G4double > ( tempLH.second , TH.second.isoAngle[ i ] ) );  
               anE_isoAng.isoAngle.push_back( angle ); 
            }
            mu = getMu( &anE_isoAng );

         } else {
            //TL.second.n != TH.second.n
            throw G4HadronicException(__FILE__, __LINE__, "A problem is found in Thermal Scattering Data! Do not yet supported");
         }
     
         //set 
         theParticleChange.SetEnergyChange( sE );
         G4double phi = CLHEP::twopi*G4UniformRand();
         G4double sint= std::sqrt ( 1 - mu*mu );
         theParticleChange.SetMomentumChange( sint*std::cos(phi), sint*std::sin(phi), mu );

      } 
      else if ( random <= ( inelastic + theXSection->GetCoherentCrossSection( dp , theElement , theMaterial ) ) / total )
      {
         // Coherent Elastic 

         G4double E = aTrack.GetKineticEnergy();

         // T_L and T_H 
         std::map < G4double , std::vector< std::pair< G4double , G4double >* >* >::iterator it; 
         std::vector<G4double> v_temp;
         v_temp.clear();
         for ( it = coherentFSs->find( ielement )->second->begin() ; it != coherentFSs->find( ielement )->second->end() ; it++ )
         {
            v_temp.push_back( it->first );
         }

//                   T_L         T_H 
         std::pair < G4double , G4double > tempLH = find_LH ( aTemp , &v_temp );
//
//
//       For T_L anEPM_TL  and T_H anEPM_TH
//
         std::vector< std::pair< G4double , G4double >* >* pvE_p_TL = NULL; 
         std::vector< std::pair< G4double , G4double >* >* pvE_p_TH = NULL; 

         if ( tempLH.first != 0.0 && tempLH.second != 0.0 ) 
         {
            pvE_p_TL = coherentFSs->find( ielement )->second->find ( tempLH.first/kelvin )->second;
            pvE_p_TH = coherentFSs->find( ielement )->second->find ( tempLH.first/kelvin )->second;
         }
         else if ( tempLH.first == 0.0 )
         {
            pvE_p_TL = coherentFSs->find( ielement )->second->find ( v_temp[ 0 ] )->second;
            pvE_p_TH = coherentFSs->find( ielement )->second->find ( v_temp[ 1 ] )->second;
            tempLH.first = tempLH.second;
            tempLH.second = v_temp[ 1 ];
         }
         else if ( tempLH.second == 0.0 )
         {
            pvE_p_TH = coherentFSs->find( ielement )->second->find ( v_temp.back() )->second;
            std::vector< G4double >::iterator itv;
            itv = v_temp.end();
            itv--;
            itv--;
            pvE_p_TL = coherentFSs->find( ielement )->second->find ( *itv )->second;
            tempLH.second = tempLH.first;
            tempLH.first = *itv;
         }
         else 
         {
            //tempLH.first == 0.0 && tempLH.second
            throw G4HadronicException(__FILE__, __LINE__, "A problem is found in Thermal Scattering Data! Unexpected temperature values in data");
         }

         std::vector< G4double > vE_T;
         std::vector< G4double > vp_T;

         G4int n1 = pvE_p_TL->size();  

         //171005 fix bug, contribution from H.N. TRAN@CEA
         for ( G4int i=0 ; i < n1 ; i++ ) 
         {
            if ( (*pvE_p_TL)[i]->first != (*pvE_p_TH)[i]->first ) throw G4HadronicException(__FILE__, __LINE__, "A problem is found in Thermal Scattering Data!");
            vE_T.push_back ( (*pvE_p_TL)[i]->first );
            vp_T.push_back ( get_linear_interpolated ( aTemp , std::pair< G4double , G4double > ( tempLH.first , (*pvE_p_TL)[i]->second ) , std::pair< G4double , G4double > ( tempLH.second , (*pvE_p_TL)[i]->second ) ) );  
         }

         G4int j = 0;  
         for ( G4int i = 1 ; i < n1 ; i++ ) 
         {
            if ( E/eV < vE_T[ i ] ) 
            {
               j = i-1;
               break;
            }
         }

         G4double rand_for_mu = G4UniformRand();

         G4int k = 0;
         for ( G4int i = 0 ; i <= j ; i++ )
         {
             G4double Pi = vp_T[ i ] / vp_T[ j ]; 
             if ( rand_for_mu < Pi )
             {
                k = i; 
                break;
             }
         }

         G4double Ei = vE_T[ k ];

         G4double mu = 1 - 2 * Ei / (E/eV) ;  
         //111102
         if ( mu < -1.0 ) mu = -1.0;
         //G4cout << "E= " << E/eV << ", Ei= " << Ei << ", mu= " << mu << G4endl;

         theParticleChange.SetEnergyChange( E );
         G4double phi = CLHEP::twopi*G4UniformRand();
         G4double sint= std::sqrt ( 1 - mu*mu );
         theParticleChange.SetMomentumChange( sint*std::cos(phi), sint*std::sin(phi), mu );

      }
      else
      {
         // InCoherent Elastic

         // T_L and T_H 
         std::map < G4double , std::vector < E_isoAng* >* >::iterator it; 
         std::vector<G4double> v_temp;
         v_temp.clear();
         for ( it = incoherentFSs->find( ielement )->second->begin() ; it != incoherentFSs->find( ielement )->second->end() ; it++ )
         {
            v_temp.push_back( it->first );
         }
              
//                   T_L         T_H 
         std::pair < G4double , G4double > tempLH = find_LH ( aTemp , &v_temp );

//
//       For T_L anEPM_TL  and T_H anEPM_TH
//

         E_isoAng anEPM_TL_E;
         E_isoAng anEPM_TH_E;

         if ( tempLH.first != 0.0 && tempLH.second != 0.0 ) {
            //Interpolate TL and TH 
            anEPM_TL_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( tempLH.first/kelvin )->second );
            anEPM_TH_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( tempLH.second/kelvin )->second );
         } else if ( tempLH.first == 0.0 ) {
            //Extrapolate T0 and T1
            anEPM_TL_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( v_temp[ 0 ] )->second );
            anEPM_TH_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( v_temp[ 1 ] )->second );
            tempLH.first = tempLH.second;
            tempLH.second = v_temp[ 1 ];
         } else if (  tempLH.second == 0.0 ) {
            //Extrapolate Tmax-1 and Tmax
            anEPM_TH_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( v_temp.back() )->second );
            std::vector< G4double >::iterator itv;
            itv = v_temp.end();
            itv--;
            itv--;
            anEPM_TL_E = create_E_isoAng_from_energy ( aTrack.GetKineticEnergy() , incoherentFSs->find( ielement )->second->find ( *itv )->second );
            tempLH.second = tempLH.first;
            tempLH.first = *itv;
         } 
        
         // E_isoAng for aTemp and aTrack.GetKineticEnergy() 
         G4double mu=1.0;
         E_isoAng anEPM_T_E;  

         if ( anEPM_TL_E.n == anEPM_TH_E.n ) 
         {
            anEPM_T_E.n = anEPM_TL_E.n; 
            for ( G4int i=0 ; i < anEPM_TL_E.n ; i++ )
            { 
               G4double angle;
               angle = get_linear_interpolated ( aTemp , std::pair< G4double , G4double > ( tempLH.first , anEPM_TL_E.isoAngle[ i ] ) , std::pair< G4double , G4double > ( tempLH.second , anEPM_TH_E.isoAngle[ i ] ) );  
               anEPM_T_E.isoAngle.push_back( angle ); 
            }
            mu = getMu ( &anEPM_T_E );

         } else {
            // anEPM_TL_E.n != anEPM_TH_E.n
            throw G4HadronicException(__FILE__, __LINE__, "A problem is found in Thermal Scattering Data! Do not yet supported");
         }

         // Set Final State
         theParticleChange.SetEnergyChange( aTrack.GetKineticEnergy() );  // No energy change in Elastic
         G4double phi = CLHEP::twopi*G4UniformRand();
         G4double sint= std::sqrt ( 1 - mu*mu );
         theParticleChange.SetMomentumChange( sint*std::cos(phi), sint*std::sin(phi), mu );

      } 
      delete dp;

      return &theParticleChange;
      
   }
   else 
   {
      // Not thermal element   
      // Neutron HP will handle
      return theHPElastic -> ApplyYourself( aTrack, aNucleus );
   }

}



G4double G4ParticleHPThermalScattering::getMu( E_isoAng* anEPM  )
{

   G4double random = G4UniformRand();
   G4double result = 0.0;  

   G4int in = int ( random * ( (*anEPM).n ) );

   if ( in != 0 )
   {
       G4double mu_l = (*anEPM).isoAngle[ in-1 ]; 
       G4double mu_h = (*anEPM).isoAngle[ in ]; 
       result = ( mu_h - mu_l ) * ( random * ( (*anEPM).n ) - in ) + mu_l; 
   }
   else 
   {
       G4double x = random * (*anEPM).n;
       //Bugzilla 1971 
       G4double ratio = 0.5;
       G4double xx = G4UniformRand();
       if ( x <= ratio ) 
       {
          G4double mu_l = -1; 
          G4double mu_h = (*anEPM).isoAngle[ 0 ]; 
          result = ( mu_h - mu_l ) * xx + mu_l; 
       }
       else
       {
          G4double mu_l = (*anEPM).isoAngle[ (*anEPM).n - 1 ]; 
          G4double mu_h = 1;
          result = ( mu_h - mu_l ) * xx + mu_l; 
       }
   }
   return result;
} 



std::pair < G4double , G4double >  G4ParticleHPThermalScattering::find_LH ( G4double x , std::vector< G4double >* aVector )
{
   G4double LL = 0.0; 
   G4double H = 0.0; 

   // v->size() == 1 --> LL=H=v(0)
   if ( aVector->size() == 1 ) {
      LL = aVector->front();
      H = aVector->front();
   } else {
   // 1) temp < v(0) -> LL=0.0 H=v(0)
   // 2) v(i-1) < temp <= v(i) -> LL=v(i-1) H=v(i)
   // 3) v(imax) < temp -> LL=v(imax) H=0.0
      for ( std::vector< G4double >::iterator 
            it = aVector->begin() ; it != aVector->end() ; it++ ) {
         if ( x <= *it ) {
            H = *it;  
            if ( it != aVector->begin() ) {
               // 2)
               it--;
               LL = *it;
            } else {
               // 1)
               LL = 0.0;
            }
            break; 
         } 
      } 
      // 3) 
      if ( H == 0.0 ) LL = aVector->back();
   }

   return std::pair < G4double , G4double > ( LL , H ); 
}



G4double G4ParticleHPThermalScattering::get_linear_interpolated ( G4double x , std::pair< G4double , G4double > Low , std::pair< G4double , G4double > High )
{ 
   G4double y=0.0;
   if ( High.first - Low.first != 0 ) {
      y = ( High.second - Low.second ) / ( High.first - Low.first ) * ( x - Low.first ) + Low.second;
   } else { 
      if ( High.second == Low.second ) {
         y = High.second;
      } else { 
         G4cout << "G4ParticleHPThermalScattering liner interpolation err!!" << G4endl; 
      }
   }
      
   return y; 
} 



E_isoAng
G4ParticleHPThermalScattering::create_E_isoAng_from_energy(G4double energy,
                                                           std::vector<E_isoAng*>* vEPM)
{
   E_isoAng anEPM_T_E;
   std::vector<E_isoAng*>::iterator iv;

   std::vector<G4double> v_e;
   v_e.clear();
   for (iv = vEPM->begin(); iv != vEPM->end(); iv++) 
     v_e.push_back( (*iv)->energy );

   std::pair<G4double, G4double> energyLH = find_LH(energy, &v_e);
   //G4cout << " " << energy/eV << " " << energyLH.first/eV  << " " << energyLH.second/eV << G4endl;

   E_isoAng* panEPM_T_EL = 0;
   E_isoAng* panEPM_T_EH = 0;

   if (energyLH.first != 0.0 && energyLH.second != 0.0) {
     for (iv = vEPM->begin(); iv != vEPM->end(); iv++) {
       if (energyLH.first == (*iv)->energy) {
         panEPM_T_EL = *iv;
         iv++;
         panEPM_T_EH = *iv;
         break;
       }
     } 
 
   } else if (energyLH.first == 0.0) {
     panEPM_T_EL = (*vEPM)[0];
     panEPM_T_EH = (*vEPM)[1];

   } else if (energyLH.second == 0.0) {
     panEPM_T_EH = (*vEPM).back();
     iv = vEPM->end();
     iv--; 
     iv--; 
     panEPM_T_EL = *iv;
   } 

   if (panEPM_T_EL != 0 && panEPM_T_EH != 0) {
     //checking isoAng has proper values or not 
     // Inelastic/FS, the first and last entries of *vEPM has all zero values.
     if ( !(check_E_isoAng(panEPM_T_EL) ) ) panEPM_T_EL = panEPM_T_EH;
     if ( !(check_E_isoAng(panEPM_T_EH) ) ) panEPM_T_EH = panEPM_T_EL;

     if (panEPM_T_EL->n == panEPM_T_EH->n) {
       anEPM_T_E.energy = energy; 
       anEPM_T_E.n = panEPM_T_EL->n; 

       for (G4int i=0; i < panEPM_T_EL->n; i++) { 
         G4double angle;
         angle = get_linear_interpolated(energy, std::pair<G4double,G4double>(energyLH.first, panEPM_T_EL->isoAngle[i] ),
                                         std::pair<G4double,G4double>(energyLH.second, panEPM_T_EH->isoAngle[i] ) );  
         anEPM_T_E.isoAngle.push_back(angle); 
       }

     } else {
       G4Exception("G4ParticleHPThermalScattering::create_E_isoAng_from_energy",
                   "NotSupported", JustWarning,
                   "G4ParticleHPThermalScattering does not support yet EL->n != EH->n."); 
     }

   } else {
     G4Exception("G4ParticleHPThermalScattering::create_E_isoAng_from_energy",
                 "HAD_THERM_000", FatalException,
                 "Pointer panEPM_T_EL or panEPM_T_EH is zero");
   }

   return anEPM_T_E;
}



G4double G4ParticleHPThermalScattering::get_secondary_energy_from_E_P_E_isoAng ( G4double random , E_P_E_isoAng* anE_P_E_isoAng )
{

   G4double secondary_energy = 0.0;

   G4int n = anE_P_E_isoAng->n;
   G4double sum_p = 0.0; // sum_p_H
   G4double sum_p_L = 0.0;

   G4double total=0.0;

/*
   delete for speed up
   for ( G4int i = 0 ; i < n-1 ; i++ )
   {
      G4double E_L = anE_P_E_isoAng->vE_isoAngle[i]->energy/eV;
      G4double E_H = anE_P_E_isoAng->vE_isoAngle[i+1]->energy/eV;
      G4double dE = E_H - E_L;
      total += ( ( anE_P_E_isoAng->prob[i] ) * dE );
   }

   if ( std::abs( total - anE_P_E_isoAng->sum_of_probXdEs ) > 1.0e-14 ) G4cout << total - anE_P_E_isoAng->sum_of_probXdEs << G4endl;
*/
   total =  anE_P_E_isoAng->sum_of_probXdEs;

   for ( G4int i = 0 ; i < n-1 ; i++ )
   {
      G4double E_L = anE_P_E_isoAng->vE_isoAngle[i]->energy/eV;
      G4double E_H = anE_P_E_isoAng->vE_isoAngle[i+1]->energy/eV;
      G4double dE = E_H - E_L;
      sum_p += ( ( anE_P_E_isoAng->prob[i] ) * dE );

      if ( random <= sum_p/total )
      {
         secondary_energy = get_linear_interpolated ( random , std::pair < G4double , G4double > ( sum_p_L/total , E_L ) , std::pair < G4double , G4double > ( sum_p/total , E_H ) );
         secondary_energy = secondary_energy*eV;  //need eV
         break;
      }
      sum_p_L = sum_p; 
   }

   return secondary_energy; 
}



std::pair< G4double , E_isoAng > G4ParticleHPThermalScattering::create_sE_and_EPM_from_pE_and_vE_P_E_isoAng ( G4double rand_for_sE ,  G4double pE , std::vector < E_P_E_isoAng* >*  vNEP_EPM )
{

         std::map< G4double , G4int > map_energy;
         map_energy.clear();
         std::vector< G4double > v_energy;
         v_energy.clear();
         std::vector< E_P_E_isoAng* >::iterator itv;
         G4int i = 0;
         for ( itv = vNEP_EPM->begin(); itv != vNEP_EPM->end(); itv++ )
         {
            v_energy.push_back( (*itv)->energy );
            map_energy.insert( std::pair < G4double , G4int > ( (*itv)->energy , i ) );
            i++;
         } 
            
         std::pair < G4double , G4double > energyLH = find_LH ( pE , &v_energy );

         E_P_E_isoAng* pE_P_E_isoAng_EL = 0; 
         E_P_E_isoAng* pE_P_E_isoAng_EH = 0; 

         if ( energyLH.first != 0.0 && energyLH.second != 0.0 ) 
         {
            pE_P_E_isoAng_EL = (*vNEP_EPM)[ map_energy.find ( energyLH.first )->second ];    
            pE_P_E_isoAng_EH = (*vNEP_EPM)[ map_energy.find ( energyLH.second )->second ];    
         }
         else if ( energyLH.first == 0.0 ) 
         {
            pE_P_E_isoAng_EL = (*vNEP_EPM)[ 0 ];    
            pE_P_E_isoAng_EH = (*vNEP_EPM)[ 1 ];    
         }
         if ( energyLH.second == 0.0 ) 
         {
            pE_P_E_isoAng_EH = (*vNEP_EPM).back();    
            itv = vNEP_EPM->end();
            itv--; 
            itv--;
            pE_P_E_isoAng_EL = *itv;    
         }


         G4double sE; 
         G4double sE_L; 
         G4double sE_H; 
         

         sE_L = get_secondary_energy_from_E_P_E_isoAng ( rand_for_sE , pE_P_E_isoAng_EL );
         sE_H = get_secondary_energy_from_E_P_E_isoAng ( rand_for_sE , pE_P_E_isoAng_EH );

         sE = get_linear_interpolated ( pE , std::pair < G4double , G4double > ( energyLH.first , sE_L ) , std::pair < G4double , G4double > ( energyLH.second , sE_H ) );  

          
         E_isoAng E_isoAng_L = create_E_isoAng_from_energy ( sE , &(pE_P_E_isoAng_EL->vE_isoAngle) );
         E_isoAng E_isoAng_H = create_E_isoAng_from_energy ( sE , &(pE_P_E_isoAng_EH->vE_isoAngle) );

         E_isoAng anE_isoAng; 
         //For defeating warning message from compiler
         anE_isoAng.n = 1;
         anE_isoAng.energy = sE; //never used 
         if ( E_isoAng_L.n == E_isoAng_H.n ) 
         {
            anE_isoAng.n =  E_isoAng_L.n; 
            for ( G4int j=0 ; j < anE_isoAng.n ; j++ )
            { 
               G4double angle;
               angle = get_linear_interpolated ( sE  , std::pair< G4double , G4double > ( sE_L , E_isoAng_L.isoAngle[ j ] ) , std::pair< G4double , G4double > ( sE_H , E_isoAng_H.isoAngle[ j ] ) );  
               anE_isoAng.isoAngle.push_back( angle ); 
            }
         }
         else
         {
            //G4cout << "Do not Suuport yet." << G4endl; 
            throw G4HadronicException(__FILE__, __LINE__, "Unexpected values!");
         }
     
   
         
   return std::pair< G4double , E_isoAng >( sE , anE_isoAng); 
}

void G4ParticleHPThermalScattering::buildPhysicsTable()
{

   //Is rebuild of physics table a necessity 
   if ( nMaterial == G4Material::GetMaterialTable()->size() && nElement == G4Element::GetElementTable()->size() ) {
      return;
   } else {
      nMaterial = G4Material::GetMaterialTable()->size(); 
      nElement = G4Element::GetElementTable()->size(); 
   }

   dic.clear();   
   std::map < G4String , G4int > co_dic;   

   //Searching Nist Materials
   static G4ThreadLocal G4MaterialTable* theMaterialTable  = 0 ; if (!theMaterialTable) theMaterialTable= G4Material::GetMaterialTable();
   size_t numberOfMaterials = G4Material::GetNumberOfMaterials();
   for ( size_t i = 0 ; i < numberOfMaterials ; i++ )
   {
      G4Material* material = (*theMaterialTable)[i];
      size_t numberOfElements = material->GetNumberOfElements();
      for ( size_t j = 0 ; j < numberOfElements ; j++ )
      {
         const G4Element* element = material->GetElement(j);
         if ( names.IsThisThermalElement ( material->GetName() , element->GetName() ) )
         {                                    
            G4int ts_ID_of_this_geometry; 
            G4String ts_ndl_name = names.GetTS_NDL_Name( material->GetName() , element->GetName() ); 
            if ( co_dic.find ( ts_ndl_name ) != co_dic.end() )
            {
               ts_ID_of_this_geometry = co_dic.find ( ts_ndl_name ) -> second;
            }
            else
            {
               ts_ID_of_this_geometry = co_dic.size();
               co_dic.insert ( std::pair< G4String , G4int >( ts_ndl_name , ts_ID_of_this_geometry ) );
            }

            //G4cout << "Neutron HP Thermal Scattering: Registering a material-element pair of " 
            //       << material->GetName() << " " << element->GetName() 
            //       << " as internal thermal scattering id of  " <<  ts_ID_of_this_geometry << "." << G4endl;

            dic.insert( std::pair < std::pair < G4Material* , const G4Element* > , G4int > ( std::pair < G4Material* , const G4Element* > ( material , element ) , ts_ID_of_this_geometry ) );
         }
      }
   }

   //Searching TS Elements 
   static G4ThreadLocal G4ElementTable* theElementTable  = 0 ; if (!theElementTable) theElementTable= G4Element::GetElementTable();
   size_t numberOfElements = G4Element::GetNumberOfElements();
   //size_t numberOfThermalElements = 0; 
   for ( size_t i = 0 ; i < numberOfElements ; i++ )
   {
      const G4Element* element = (*theElementTable)[i];
      if ( names.IsThisThermalElement ( element->GetName() ) )
      {
         if ( names.IsThisThermalElement ( element->GetName() ) )
         {                                    
            G4int ts_ID_of_this_geometry; 
            G4String ts_ndl_name = names.GetTS_NDL_Name( element->GetName() ); 
            if ( co_dic.find ( ts_ndl_name ) != co_dic.end() )
            {
               ts_ID_of_this_geometry = co_dic.find ( ts_ndl_name ) -> second;
            }
            else
            {
               ts_ID_of_this_geometry = co_dic.size();
               co_dic.insert ( std::pair< G4String , G4int >( ts_ndl_name , ts_ID_of_this_geometry ) );
            }

            //G4cout << "Neutron HP Thermal Scattering: Registering an element of " 
            //       << material->GetName() << " " << element->GetName() 
            //       << " as internal thermal scattering id of  " <<  ts_ID_of_this_geometry << "." << G4endl;

            dic.insert( std::pair < std::pair < const G4Material* , const G4Element* > , G4int > ( std::pair < const G4Material* , const G4Element* > ( (G4Material*)NULL , element ) ,  ts_ID_of_this_geometry ) );
         }
      }
   }

   G4cout << G4endl;
   G4cout << "Neutron HP Thermal Scattering: Following material-element pairs or elements are registered." << G4endl;
   for ( std::map < std::pair < const G4Material* , const G4Element* > , G4int >::iterator it = dic.begin() ; it != dic.end() ; it++ )   
   {
      if ( it->first.first != NULL ) 
      {
         G4cout << "Material " << it->first.first->GetName() << " - Element " << it->first.second->GetName() << ",  internal thermal scattering id " << it->second << G4endl;
      }
      else
      {
         G4cout << "Element " << it->first.second->GetName() << ",  internal thermal scattering id " << it->second << G4endl;
      }
   }
   G4cout << G4endl;

   // Read Cross Section Data files
   
   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();
   coherentFSs = hpmanager->GetThermalScatteringCoherentFinalStates();
   incoherentFSs = hpmanager->GetThermalScatteringIncoherentFinalStates();
   inelasticFSs = hpmanager->GetThermalScatteringInelasticFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      clearCurrentFSData();

      if ( coherentFSs == NULL ) coherentFSs = new std::map < G4int , std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* >;
      if ( incoherentFSs == NULL ) incoherentFSs = new std::map < G4int , std::map < G4double , std::vector < E_isoAng* >* >* >;
      if ( inelasticFSs == NULL ) inelasticFSs = new std::map < G4int , std::map < G4double , std::vector < E_P_E_isoAng* >* >* >;

       G4String dirName;
       if ( !std::getenv( "G4NEUTRONHPDATA" ) ) 
          throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
       dirName = std::getenv( "G4NEUTRONHPDATA" );

   //G4String name;

   for ( std::map < G4String , G4int >::iterator it = co_dic.begin() ; it != co_dic.end() ; it++ )  
   {
      G4String tsndlName = it->first;
      G4int ts_ID = it->second;

      // Coherent
      G4String fsName = "/ThermalScattering/Coherent/FS/";
      G4String fileName = dirName + fsName + tsndlName;
      coherentFSs->insert ( std::pair < G4int , std::map < G4double , std::vector < std::pair< G4double , G4double >* >* >* > ( ts_ID , readACoherentFSDATA( fileName ) ) ); 

      // incoherent elastic 
      fsName = "/ThermalScattering/Incoherent/FS/";
      fileName = dirName + fsName + tsndlName;
      incoherentFSs->insert ( std::pair < G4int , std::map < G4double , std::vector < E_isoAng* >* >* > ( ts_ID , readAnIncoherentFSDATA( fileName ) ) ); 

      // inelastic 
      fsName = "/ThermalScattering/Inelastic/FS/";
      fileName = dirName + fsName + tsndlName;
      inelasticFSs->insert ( std::pair < G4int , std::map < G4double , std::vector < E_P_E_isoAng* >* >* > ( ts_ID , readAnInelasticFSDATA( fileName ) ) ); 
   } 

      hpmanager->RegisterThermalScatteringCoherentFinalStates( coherentFSs );
      hpmanager->RegisterThermalScatteringIncoherentFinalStates( incoherentFSs );
      hpmanager->RegisterThermalScatteringInelasticFinalStates( inelasticFSs );
   }

   theXSection->BuildPhysicsTable( *(G4Neutron::Neutron()) );
}
 

G4int G4ParticleHPThermalScattering::getTS_ID ( const G4Material* material , const G4Element* element )
{
   G4int result = -1;
   if ( dic.find( std::pair < const G4Material* , const G4Element* > ( material , element ) ) != dic.end() ) 
      result = dic.find( std::pair < const G4Material* , const G4Element* > ( material , element ) )->second; 
   return result; 
}

const std::pair<G4double, G4double> G4ParticleHPThermalScattering::GetFatalEnergyCheckLevels() const
{
   //return std::pair<G4double, G4double>(10*perCent,10*GeV);
   return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}

void G4ParticleHPThermalScattering::AddUserThermalScatteringFile( G4String nameG4Element , G4String filename)
{
   names.AddThermalElement( nameG4Element , filename );
   theXSection->AddUserThermalScatteringFile( nameG4Element , filename );
   buildPhysicsTable();
}


G4bool G4ParticleHPThermalScattering::check_E_isoAng( E_isoAng* anE_IsoAng )
{
   G4bool result=false;

   G4int n = anE_IsoAng->n;
   G4double sum=0.0;
   for ( G4int i = 0 ; i < n ; i++ ) {
      sum += anE_IsoAng->isoAngle[ i ];
   }
   if ( sum != 0.0 ) result = true;

   return result;
}

void G4ParticleHPThermalScattering::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on thermal scattering data in\n"
           << "evaluated nuclear data libraries for neutrons below 5eV\n"
           << "on specific materials\n";
}
