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
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPCaptureFS.hh"
#include "G4ParticleHPDeExGammas.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Threading.hh"

  G4ParticleHPCapture::G4ParticleHPCapture()
   :G4HadronicInteraction("NeutronHPCapture")
  ,theCapture(NULL)
  ,numEle(0)
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
/*
//    G4cout << "Capture : start of construction!!!!!!!!"<<G4endl;
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Capture";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
//    G4cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
//    G4cout <<"Disname="<<dirName<<" numEle="<<numEle<<G4endl;
    //theCapture = new G4ParticleHPChannel[numEle];
//    G4cout <<"G4ParticleHPChannel constructed"<<G4endl;
    G4ParticleHPCaptureFS * theFS = new G4ParticleHPCaptureFS;
    //for (G4int i=0; i<numEle; i++)
    //{
//  //    G4cout << "initializing theCapture "<<i<<" "<< numEle<<G4endl;
    //  theCapture[i].Init((*(G4Element::GetElementTable()))[i], dirName);
    //  theCapture[i].Register(theFS);
    //}
    for ( G4int i = 0 ; i < numEle ; i++ ) 
    {
       theCapture.push_back( new G4ParticleHPChannel );
       (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
       (*theCapture[i]).Register(theFS);
    }
    delete theFS;
//    G4cout << "-------------------------------------------------"<<G4endl;
//    G4cout << "Leaving G4ParticleHPCapture::G4ParticleHPCapture"<<G4endl;
*/
  }
  
  G4ParticleHPCapture::~G4ParticleHPCapture()
  {
    //delete [] theCapture;
    //vector is shared, only master deletes
    if ( ! G4Threading::IsWorkerThread() ) {
        if ( theCapture != NULL ) {
            for ( std::vector<G4ParticleHPChannel*>::iterator
                ite = theCapture->begin() ; ite != theCapture->end() ; ite++ ) {
                delete *ite;
            }
            theCapture->clear();
        }
    }
  }
  
  #include "G4ParticleHPThermalBoost.hh"
  G4HadFinalState * G4ParticleHPCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
  {

    //if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();

    G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();
    if(getenv("NeutronHPCapture")) G4cout <<" ####### G4ParticleHPCapture called"<<G4endl;
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      G4double* xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4ParticleHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        //xSec[i] = theCapture[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
        xSec[i] = ((*theCapture)[index])->GetXsec(aThermalE.GetThermalEnergy(aTrack,
  		                                                     theMaterial->GetElement(i),
  								     theMaterial->GetTemperature()));
        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; i++)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        //if(random<=running/sum) break;
        if( sum == 0 || random <= running/sum ) break;
      }
      if(i==n) i=std::max(0, n-1);
      delete [] xSec;
    }

    //return theCapture[index].ApplyYourself(aTrack);
    //G4HadFinalState* result = theCapture[index].ApplyYourself(aTrack);
    G4HadFinalState* result = ((*theCapture)[index])->ApplyYourself(aTrack);

    //Overwrite target parameters
    aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; j++ ) { 
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
    }
    //G4cout << "Target Material of this reaction is " << theMaterial->GetName() << G4endl;
    //G4cout << "Target Element of this reaction is " << target_element->GetName() << G4endl;
    //G4cout << "Target Isotope of this reaction is " << target_isotope->GetName() << G4endl;
    aNucleus.SetIsotope( target_isotope );

    G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
    return result; 
  }

const std::pair<G4double, G4double> G4ParticleHPCapture::GetFatalEnergyCheckLevels() const
{
   //return std::pair<G4double, G4double>(10*perCent,10*GeV);
   return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}

/*
void G4ParticleHPCapture::addChannelForNewElement()
{
   G4ParticleHPCaptureFS* theFS = new G4ParticleHPCaptureFS;
   for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
   {
      G4cout << "G4ParticleHPCapture Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;
      theCapture.push_back( new G4ParticleHPChannel );
      (*theCapture[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
      (*theCapture[i]).Register(theFS);
   }
   delete theFS;
   numEle = (G4int)G4Element::GetNumberOfElements();
}
*/

G4int G4ParticleHPCapture::GetVerboseLevel() const
{
   return G4ParticleHPManager::GetInstance()->GetVerboseLevel();
}
void G4ParticleHPCapture::SetVerboseLevel( G4int newValue ) 
{
   G4ParticleHPManager::GetInstance()->SetVerboseLevel(newValue);
}

void G4ParticleHPCapture::BuildPhysicsTable(const G4ParticleDefinition&)
{

   G4ParticleHPManager* hpmanager = G4ParticleHPManager::GetInstance();

   theCapture = hpmanager->GetCaptureFinalStates();

   if ( G4Threading::IsMasterThread() ) {

      if ( theCapture == NULL ) theCapture = new std::vector<G4ParticleHPChannel*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theCapture->size() == G4Element::GetNumberOfElements() ) {
         numEle = G4Element::GetNumberOfElements();
         return;
      }

      if ( !getenv("G4NEUTRONHPDATA") ) 
          throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = getenv("G4NEUTRONHPDATA");
      G4String tString = "/Capture";
      dirName = dirName + tString;

      G4ParticleHPCaptureFS * theFS = new G4ParticleHPCaptureFS;
      for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
      {
         theCapture->push_back( new G4ParticleHPChannel );
         ((*theCapture)[i])->Init((*(G4Element::GetElementTable()))[i], dirName);
         ((*theCapture)[i])->Register(theFS);
      }
      delete theFS;
      hpmanager->RegisterCaptureFinalStates( theCapture );
   }
   numEle = G4Element::GetNumberOfElements();
}

void G4ParticleHPCapture::ModelDescription(std::ostream& outFile) const
{
   outFile << "High Precision model based on Evaluated Nuclear Data Files (ENDF) for radiative capture reaction of neutrons below 20MeV\n";
}
