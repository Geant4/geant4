//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPorLCapture.hh"
#include "G4NeutronHPCaptureFS.hh"

G4NeutronHPorLCapture::G4NeutronHPorLCapture()
{
   G4NeutronHPCaptureFS * theFS = new G4NeutronHPCaptureFS;
   if(!getenv("NeutronHPCrossSections")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
   dirName = getenv("NeutronHPCrossSections");
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
      try { while(!theCapture[i].Register(theFS)); }
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
