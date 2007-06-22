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
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureFS.hh"
#include "G4NeutronHPDeExGammas.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

  G4NeutronHPCapture::G4NeutronHPCapture()
   :G4HadronicInteraction("NeutronHPCapture")
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
//    G4cout << "Capture : start of construction!!!!!!!!"<<G4endl;
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Capture/";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
//    G4cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"<<G4endl;
//    G4cout <<"Disname="<<dirName<<" numEle="<<numEle<<G4endl;
    theCapture = new G4NeutronHPChannel[numEle];
//    G4cout <<"G4NeutronHPChannel constructed"<<G4endl;
    G4NeutronHPCaptureFS * theFS = new G4NeutronHPCaptureFS;
    for (G4int i=0; i<numEle; i++)
    {
//      G4cout << "initializing theCapture "<<i<<" "<< numEle<<G4endl;
      theCapture[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      theCapture[i].Register(theFS);
    }
    delete theFS;
//    G4cout << "-------------------------------------------------"<<G4endl;
//    G4cout << "Leaving G4NeutronHPCapture::G4NeutronHPCapture"<<G4endl;
  }
  
  G4NeutronHPCapture::~G4NeutronHPCapture()
  {
    delete [] theCapture;
//    G4cout << "Leaving G4NeutronHPCapture::~G4NeutronHPCapture"<<G4endl;
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  G4HadFinalState * G4NeutronHPCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& )
  {
    if(getenv("NeutronHPCapture")) G4cout <<" ####### G4NeutronHPCapture called"<<G4endl;
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    if(n!=1)
    {
      xSec = new G4double[n];
      G4double sum=0;
      G4int i;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        xSec[i] = theCapture[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
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
    return theCapture[index].ApplyYourself(aTrack);
  }
