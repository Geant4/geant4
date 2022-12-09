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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPFFFissionFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"

G4ParticleHPFFFissionFS::~G4ParticleHPFFFissionFS()
{
    std::map<G4int,std::map<G4double,std::map<G4int,G4double >* >* >::iterator it = FissionProductYieldData.begin();
    while ( it != FissionProductYieldData.end() ) { // Loop checking, 11.05.2015, T. Koi
        std::map<G4double,std::map<G4int,G4double>* >* firstLevel = it->second;
        if ( firstLevel ) {
            std::map<G4double,std::map<G4int,G4double>*>::iterator it2 = firstLevel->begin();
            while ( it2 != firstLevel->end() ) { // Loop checking, 11.05.2015, T. Koi
                delete it2->second;
                it2->second = 0;
                firstLevel->erase(it2);
                it2=firstLevel->begin();
            }
        }
        delete firstLevel;
        it->second = 0;
        FissionProductYieldData.erase(it);
        it = FissionProductYieldData.begin();
    }
    
    std::map< G4int , std::map< G4double , G4int >* >::iterator ii = mMTInterpolation.begin();
    while ( ii != mMTInterpolation.end() ) { // Loop checking, 11.05.2015, T. Koi
        delete ii->second;
        mMTInterpolation.erase(ii);
        ii = mMTInterpolation.begin();
    }
}

void G4ParticleHPFFFissionFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & , G4ParticleDefinition*)
{
   //G4cout << "G4ParticleHPFFFissionFS::Init" << G4endl;
   G4String aString = "FF";

   G4String tString = dirName;
   G4bool dbool;
   G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, tString, aString , dbool);
   G4String filename = aFile.GetName();
   theBaseA = aFile.GetA();
   theBaseZ = aFile.GetZ();

//3456
   if ( !dbool || ( Z < 2.5 && ( std::abs(theBaseZ-Z)>0.0001 || std::abs(theBaseA-A)>0.0001) ) )
   {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return; // no data for exactly this isotope.
   }
   //std::ifstream theData(filename, std::ios::in);
   std::istringstream theData(std::ios::in);
   G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);
   G4double dummy;
   if ( !theData )
   {
      //theData.close();
      hasFSData = false;
      hasXsec = false;
      hasAnyData = false;
      return; // no data for this FS for this isotope
   }


   hasFSData = true; 
      //          MT              Energy            FPS    Yield
      //std::map< int , std::map< double , std::map< int , double >* >* > FisionProductYieldData; 
   while ( theData.good() ) // Loop checking, 11.05.2015, T. Koi
   {
      G4int iMT, iMF;
      G4int imax;
      //Reading the data
      //         MT       MF       AWR
      theData >> iMT >> iMF >> dummy;
      //         nBlock 
      theData >> imax;
      //if ( !theData.good() ) continue;
      //        Ei                   FPS     Yield
      std::map< G4double , std::map< G4int , G4double >* >* mEnergyFSPData = new std::map< G4double , std::map< G4int , G4double >* >;

      std::map< G4double , G4int >* mInterporation = new std::map< G4double , G4int >;
      for ( G4int i = 0 ; i <= imax ; i++ ) 
      {
       
         G4double YY=0.0; 
         G4double Ei;
         G4int jmax;
         G4int ip;
         //        energy of incidence neutron
         theData >> Ei;
         //        Number of data set followings 
         theData >> jmax;
         //         interpolation scheme
         theData >> ip;
         mInterporation->insert( std::pair<G4double,G4int>(Ei*eV,ip) );
         //         nNumber  nIP
         std::map<G4int,G4double>* mFSPYieldData = new std::map<G4int,G4double>;
         for ( G4int j = 0 ; j < jmax ; j++ ) 
         {
            G4int FSP;
            G4int mFSP;
            G4double Y;
            theData >> FSP >> mFSP >> Y;
            G4int k = FSP*100+mFSP;
            YY = YY + Y;
            //if ( iMT == 454 )G4cout << iMT << " " << i << " " << j << " " <<  k << " " << Y << " " << YY << G4endl;
            mFSPYieldData->insert( std::pair<G4int,G4double>( k , YY ) );
         }
         mEnergyFSPData->insert( std::pair<G4double,std::map<G4int,G4double>*>(Ei*eV,mFSPYieldData) ); 
      }

      FissionProductYieldData.insert( std::pair< G4int , std::map< G4double , std::map< G4int , G4double >* >* > (iMT,mEnergyFSPData));
      mMTInterpolation.insert( std::pair<G4int,std::map<G4double,G4int>*> (iMT,mInterporation) ); 
   } 
   //theData.close();
}
  
G4DynamicParticleVector * G4ParticleHPFFFissionFS::ApplyYourself(G4int nNeutrons)
{  
   G4DynamicParticleVector * aResult;
//    G4cout <<"G4ParticleHPFFFissionFS::ApplyYourself +"<<G4endl;
   aResult = G4ParticleHPFissionBaseFS::ApplyYourself(nNeutrons);    
   return aResult;
}

void G4ParticleHPFFFissionFS::GetAFissionFragment( G4double energy , G4int& fragZ , G4int& fragA , G4int& fragM )
{
   //G4cout << "G4ParticleHPFFFissionFS::GetAFissionFragment " << G4endl;

   G4double rand =G4UniformRand();
   //G4cout << rand << G4endl;

   std::map< G4double , std::map< G4int , G4double >* >* mEnergyFSPData = FissionProductYieldData.find( 454 )->second;
    
   //It is not clear that the treatment of the scheme 2 on two-dimensional interpolation. 
   //So, here just use the closest energy point array of yield data. 
   //TK120531
   G4double key_energy = DBL_MAX;
   if ( mEnergyFSPData->size() == 1 )
   {
      key_energy = mEnergyFSPData->cbegin()->first;
   }
   else
   {
      //Find closest energy point
      G4double Dmin=DBL_MAX; 
      for ( auto it = mEnergyFSPData->cbegin(); it != mEnergyFSPData->cend(); ++it )
      {
         G4double e = (it->first);
         G4double d = std::fabs ( energy - e ); 
         if ( d < Dmin ) 
         {
            Dmin = d;
            key_energy = e;
         }
      }
   }

   std::map<G4int,G4double>* mFSPYieldData = (*mEnergyFSPData)[key_energy];

   G4int ifrag=0;
   G4double ceilling = mFSPYieldData->rbegin()->second; // Because of numerical accuracy, this is not always 2
   for ( auto it = mFSPYieldData->cbegin(); it != mFSPYieldData->cend(); ++it )
   {
      //if ( ( rand - it->second/ceilling ) < 1.0e-6 ) std::cout << rand - it->second/ceilling << std::endl;
      if ( rand <= it->second/ceilling ) 
      {  
         //G4cout << it->first << " " << it->second/ceilling << G4endl;
         ifrag = it->first;
         break;
      }
   }

   fragZ = ifrag/100000;
   fragA = (ifrag%100000)/100;
   fragM = (ifrag%100);

   //G4cout << fragZ << " " << fragA << " " << fragM << G4endl;
}
