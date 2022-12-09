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
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 080319 Compilation warnings - gcc-4.3.0 fix by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPElasticFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4Threading.hh"
#include "G4ParticleHPThermalBoost.hh"


G4ParticleHPElastic::G4ParticleHPElastic() 
  : G4HadronicInteraction("NeutronHPElastic"), theElastic(nullptr), numEle(0)
{
   overrideSuspension = false;
   SetMinEnergy(0.*eV);
   SetMaxEnergy(20.*MeV);
}

  
G4ParticleHPElastic::~G4ParticleHPElastic()
{
   //the vectror is shared among threads, only master deletes
   if ( ! G4Threading::IsWorkerThread() ) {
      if ( theElastic != nullptr ) {
         for ( auto it=theElastic->cbegin(); it!=theElastic->cend(); ++it ) {
            delete *it;
         }
         theElastic->clear();
      }
   }
}


G4HadFinalState * G4ParticleHPElastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus)
{
   return this->ApplyYourself(aTrack, aNucleus, 0);
}
  

//--------------------------------------------------------
// New method added by L. Thulliez (CEA-Saclay) 2021/05/04
//--------------------------------------------------------
G4HadFinalState * G4ParticleHPElastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus, G4bool isFromTSL)
{
   G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
   const G4Material * theMaterial = aTrack.GetMaterial();
   G4int n = (G4int)theMaterial->GetNumberOfElements();
   std::size_t index = theMaterial->GetElement(0)->GetIndex();
 
   if ( ! isFromTSL ) {
      if ( n != 1 ) {
         G4int i;
         G4double* xSec = new G4double[n];
         G4double sum=0;
         const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
         G4double rWeight;
         G4ParticleHPThermalBoost aThermalE;
         for ( i = 0; i < n; ++i ) {
            index = theMaterial->GetElement(i)->GetIndex();
            rWeight = NumAtomsPerVolume[i];
            xSec[i] = ((*theElastic)[index])->GetXsec(aThermalE.GetThermalEnergy(aTrack,
                                                                                 theMaterial->GetElement(i),
                                                                                 theMaterial->GetTemperature()));
            xSec[i] *= rWeight;
            sum+=xSec[i];
         }
         G4double random = G4UniformRand();
         G4double running = 0;
         for ( i = 0; i < n; ++i ) {
           running += xSec[i];
           index = theMaterial->GetElement(i)->GetIndex();
           if ( sum == 0 || random <= running/sum ) break;
         }
         delete [] xSec;
      }
   } else {
      G4int i;
      if ( n != 1 ) {
         for ( i = 0; i < n; ++i ) {
            if ( aNucleus.GetZ_asInt() == (G4int)(theMaterial->GetElement(i)->GetZ()) ) {
               index = theMaterial->GetElement(i)->GetIndex();
            }
         }
      }
   }
 	
   G4HadFinalState* finalState = ((*theElastic)[index])->ApplyYourself(aTrack);
   if (overrideSuspension) finalState->SetStatusChange(isAlive);
 
   // Overwrite target parameters
   aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
   const G4Element* target_element = (*G4Element::GetElementTable())[index];
   const G4Isotope* target_isotope=nullptr;
   G4int iele = (G4int)target_element->GetNumberOfIsotopes();
   for ( G4int j = 0 ; j != iele ; ++j ) { 
      target_isotope=target_element->GetIsotope( j );
      if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
   }
   aNucleus.SetIsotope( target_isotope );
   
   G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
   return finalState; 
}


const std::pair<G4double, G4double> G4ParticleHPElastic::GetFatalEnergyCheckLevels() const
{
   // max energy non-conservation is mass of heavy nucleus
   return std::pair<G4double, G4double>(10.0*perCent, 350.0*CLHEP::GeV);
}


G4int G4ParticleHPElastic::GetVerboseLevel() const 
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}


void G4ParticleHPElastic::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}


void G4ParticleHPElastic::BuildPhysicsTable(const G4ParticleDefinition&)
{

   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   theElastic = hpmanager->GetElasticFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      if ( theElastic == nullptr ) theElastic = new std::vector<G4ParticleHPChannel*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theElastic->size() == G4Element::GetNumberOfElements() ) {
         numEle = (G4int)G4Element::GetNumberOfElements();
         return;
      }

      G4ParticleHPElasticFS * theFS = new G4ParticleHPElasticFS;
      if(!G4FindDataDir("G4NEUTRONHPDATA"))
         throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = G4FindDataDir("G4NEUTRONHPDATA");
      G4String tString = "/Elastic";
      dirName = dirName + tString;
      for ( G4int i = numEle; i < (G4int)G4Element::GetNumberOfElements(); ++i ) {
         theElastic->push_back( new G4ParticleHPChannel );
         ((*theElastic)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
         //while(!((*theElastic)[i])->Register(theFS)) ;
         ((*theElastic)[i])->Register(theFS) ;
      }
      delete theFS;
      hpmanager->RegisterElasticFinalStates( theElastic );

   }
   numEle = (G4int)G4Element::GetNumberOfElements();
}


void G4ParticleHPElastic::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF) for inelastic reaction of neutrons below 20MeV\n";
}
