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
//
// 080319 Compilation warnings - gcc-4.3.0 fix by T. Koi
//

// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPorLCapture.hh"
#include "G4SystemOfUnits.hh"
#include "G4NeutronHPCaptureFS.hh"

G4NeutronHPorLCapture::G4NeutronHPorLCapture()
  :G4HadronicInteraction("NeutronHPorLCapture")
{
   G4NeutronHPCaptureFS * theFS = new G4NeutronHPCaptureFS;
   if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
   dirName = getenv("G4NEUTRONHPDATA");
   G4String tString = "/Capture/";
   dirName = dirName + tString;
//    G4cout <<"G4NeutronHPorLCapture::G4NeutronHPorLCapture testit "<<dirName<<G4endl;
   numEle = G4Element::GetNumberOfElements();
   theCapture = new G4NeutronHPChannel[numEle];
   unavailable_elements.clear();
   for (G4int i=0; i<numEle; i++)
   {
      theCapture[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      //G4cout << (*(G4Element::GetElementTable()))[i] -> GetName()  << G4endl;
      //while(!theCapture[i].Register(theFS));
      try { while(!theCapture[i].Register(theFS)) ; }
      catch ( G4HadronicException )
      {
          unavailable_elements.insert ( (*(G4Element::GetElementTable()))[i]->GetName() ); 
      }
   }
   delete theFS;
   SetMinEnergy(0.*eV);
   SetMaxEnergy(20.*MeV);
   if ( unavailable_elements.size() > 0 )
   {
      std::set< G4String>::iterator it;
      G4cout << "HP Capture data are not available for thess  elements "<< G4endl;
      for ( it = unavailable_elements.begin() ; it != unavailable_elements.end() ; it++ )
         G4cout << *it << G4endl;
      G4cout << "Low Energy Parameterization Models will be used."<< G4endl;
   }

   createXSectionDataSet();
}
  
G4NeutronHPorLCapture::~G4NeutronHPorLCapture()
{
   delete [] theCapture;
   delete theDataSet; 
}
  
#include "G4NeutronHPThermalBoost.hh"
  
G4HadFinalState * G4NeutronHPorLCapture::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& )
{
   const G4Material * theMaterial = aTrack.GetMaterial();
   G4int n = theMaterial->GetNumberOfElements();
   G4int index = theMaterial->GetElement(0)->GetIndex();
   if(n!=1)
   {
      G4int i;
      xSec = new G4double[n];
      G4double sum=0;
      const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
      G4double rWeight;    
      G4NeutronHPThermalBoost aThermalE;
      for (i=0; i<n; i++)
      {
        index = theMaterial->GetElement(i)->GetIndex();
        rWeight = NumAtomsPerVolume[i];
        G4double x = aThermalE.GetThermalEnergy(aTrack, theMaterial->GetElement(i), theMaterial->GetTemperature());

        //xSec[i] = theCapture[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
        //		                                                     theMaterial->GetElement(i),
        //								     theMaterial->GetTemperature()));
        xSec[i] = theCapture[index].GetXsec(x);

        xSec[i] *= rWeight;
        sum+=xSec[i];
      }
      G4double random = G4UniformRand();
      G4double running = 0;
      for (i=0; i<n; i++)
      {
        running += xSec[i];
        index = theMaterial->GetElement(i)->GetIndex();
        if(random<=running/sum) break;
      }
      delete [] xSec;
      // it is element-wise initialised.
    }
    return theCapture[index].ApplyYourself(aTrack); 
}



G4bool G4NeutronHPorLCapture::IsThisElementOK( G4String name )
{
   if ( unavailable_elements.find( name ) == unavailable_elements.end() ) 
      return TRUE;
   else 
      return FALSE; 
}



void G4NeutronHPorLCapture::createXSectionDataSet()
{
   theDataSet = new G4NeutronHPorLCaptureData ( theCapture , &unavailable_elements );
}
const std::pair<G4double, G4double> G4NeutronHPorLCapture::GetFatalEnergyCheckLevels() const
{
   //return std::pair<G4double, G4double>(10*perCent,10*GeV);
   return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}
