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
// 08-08-06 delete unnecessary and harmed declaration; Bug Report[857]
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPFission.hh"
#include "G4ParticleHPFissionFS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPThermalBoost.hh"
#include "G4Threading.hh"

  G4ParticleHPFission::G4ParticleHPFission()
    : G4HadronicInteraction("NeutronHPFission"), theFission(nullptr), numEle(0)
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
  }
  
  G4ParticleHPFission::~G4ParticleHPFission()
  {
    // Vector is shared, only master deletes it
    // delete [] theFission;
    if ( ! G4Threading::IsMasterThread() ) {
        if ( theFission != nullptr ) {
            for ( auto it=theFission->cbegin(); it!=theFission->cend(); ++it ) {
                delete *it;
            }
            theFission->clear();
        }
    }
  }
  
  G4HadFinalState * G4ParticleHPFission::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
  {

    G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = (G4int)theMaterial->GetNumberOfElements();
    std::size_t index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      G4double* xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4ParticleHPThermalBoost aThermalE;
      for (i=0; i<n; ++i)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        xSec[i] = ((*theFission)[index])
                  ->GetXsec(aThermalE.GetThermalEnergy(aTrack,
                                                theMaterial->GetElement(i),
                                                theMaterial->GetTemperature()));
        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; ++i)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        //if(random<=running/sum) break;
        if( sum == 0 ||  random <= running/sum ) break;
      }
      delete [] xSec;
    }
    //return theFission[index].ApplyYourself(aTrack);                 //-2:Marker for Fission
    G4HadFinalState* result = ((*theFission)[index])->ApplyYourself(aTrack,-2);

    //Overwrite target parameters
    aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = (G4int)target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; ++j ) { 
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
    }
    aNucleus.SetIsotope( target_isotope );

    G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
    return result; 
  }

const std::pair<G4double, G4double> G4ParticleHPFission::GetFatalEnergyCheckLevels() const
{
   // max energy non-conservation is mass of heavy nucleus
   return std::pair<G4double, G4double>(10.0*perCent, 350.0*CLHEP::GeV);
}

G4int G4ParticleHPFission::GetVerboseLevel() const
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}

void G4ParticleHPFission::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void G4ParticleHPFission::BuildPhysicsTable(const G4ParticleDefinition&)
{

   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   theFission = hpmanager->GetFissionFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      if ( theFission == nullptr ) theFission = new std::vector<G4ParticleHPChannel*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theFission->size() == G4Element::GetNumberOfElements() ) {
         numEle = (G4int)G4Element::GetNumberOfElements();
         return;
      }

      if ( !G4FindDataDir("G4NEUTRONHPDATA") )
          throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = G4FindDataDir("G4NEUTRONHPDATA");
      G4String tString = "/Fission";
      dirName = dirName + tString;

      for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; ++i ) {
         theFission->push_back( new G4ParticleHPChannel );
         if ((*(G4Element::GetElementTable()))[i]->GetZ()>87) { //TK modified for ENDF-VII
            ((*theFission)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
            ((*theFission)[i])->Register( new G4ParticleHPFissionFS );
         }
      }
      hpmanager->RegisterFissionFinalStates( theFission );
   }
   numEle = (G4int)G4Element::GetNumberOfElements();
}

void G4ParticleHPFission::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF)\n"
           << "for induced fission reaction of neutrons below 20MeV\n";
}
