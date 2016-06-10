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
// this code implementation is the intellectual property of
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard (and T. Koi)
// 081203 limit maximum trial for creating final states add protection for 1H isotope case by T. Koi
//
#include "G4NeutronHPInelastic.hh"
#include "G4SystemOfUnits.hh"

#include "G4NeutronHPManager.hh"

  G4NeutronHPInelastic::G4NeutronHPInelastic()
    :G4HadronicInteraction("NeutronHPInelastic")
  ,theInelastic(NULL)
  ,numEle(0)
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
/*

    G4int istatus; 
#if defined WIN32-VC
    istatus = system("echo %G4NEUTRONHPDATA%");
#else
    istatus = system("echo $G4NEUTRONHPDATA");
#endif
    if ( istatus < 0 )
    {
      G4cout << "Warning! system(\"echo $G4NEUTRONHPDATA\") returns error value at G4NeutronHPInelastic" << G4endl;
    } 

//    G4cout << " entering G4NeutronHPInelastic constructor"<<G4endl;
    if(!getenv("G4NEUTRONHPDATA")) 
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
    dirName = getenv("G4NEUTRONHPDATA");
    G4String tString = "/Inelastic";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
*/
/*
    theInelastic = new G4NeutronHPChannelList[numEle];
    for (G4int i=0; i<numEle; i++)
    { 
      theInelastic[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      G4int itry = 0;
      do
      {
	theInelastic[i].Register(&theNFS, "F01"); // has
	theInelastic[i].Register(&theNXFS, "F02");
	theInelastic[i].Register(&the2NDFS, "F03");
 	theInelastic[i].Register(&the2NFS, "F04"); // has, E Done
 	theInelastic[i].Register(&the3NFS, "F05"); // has, E Done
  	theInelastic[i].Register(&theNAFS, "F06");
	theInelastic[i].Register(&theN3AFS, "F07");
	theInelastic[i].Register(&the2NAFS, "F08");
	theInelastic[i].Register(&the3NAFS, "F09");
	theInelastic[i].Register(&theNPFS, "F10");
	theInelastic[i].Register(&theN2AFS, "F11");
	theInelastic[i].Register(&the2N2AFS, "F12");
	theInelastic[i].Register(&theNDFS, "F13");
	theInelastic[i].Register(&theNTFS, "F14");
	theInelastic[i].Register(&theNHe3FS, "F15");
	theInelastic[i].Register(&theND2AFS, "F16");
	theInelastic[i].Register(&theNT2AFS, "F17");
	theInelastic[i].Register(&the4NFS, "F18"); // has, E Done
	theInelastic[i].Register(&the2NPFS, "F19");
	theInelastic[i].Register(&the3NPFS, "F20");
	theInelastic[i].Register(&theN2PFS, "F21");
	theInelastic[i].Register(&theNPAFS, "F22");
  	theInelastic[i].Register(&thePFS, "F23");
	theInelastic[i].Register(&theDFS, "F24");
	theInelastic[i].Register(&theTFS, "F25");
	theInelastic[i].Register(&theHe3FS, "F26");
	theInelastic[i].Register(&theAFS, "F27");
	theInelastic[i].Register(&the2AFS, "F28");
	theInelastic[i].Register(&the3AFS, "F29");
	theInelastic[i].Register(&the2PFS, "F30");
	theInelastic[i].Register(&thePAFS, "F31");
	theInelastic[i].Register(&theD2AFS, "F32");
	theInelastic[i].Register(&theT2AFS, "F33");
	theInelastic[i].Register(&thePDFS, "F34");
	theInelastic[i].Register(&thePTFS, "F35");
	theInelastic[i].Register(&theDAFS, "F36");
	theInelastic[i].RestartRegistration();
        itry++;
      }
      //while(!theInelastic[i].HasDataInAnyFinalState());
      while( !theInelastic[i].HasDataInAnyFinalState() && itry < 6 );
                                                              // 6 is corresponding to the value(5) of G4NeutronHPChannel. TK  

      if ( itry == 6 ) 
      {
         // No Final State at all.
         G4bool exceptional = false;
         if ( (*(G4Element::GetElementTable()))[i]->GetNumberOfIsotopes() == 1 )
         {
            if ( (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetZ() == 1 && (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetN() == 1 ) exceptional = true;  //1H
         } 
         if ( !exceptional ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this element");
      }
    }
*/
/*

    for (G4int i=0; i<numEle; i++)
    { 
      theInelastic.push_back( new G4NeutronHPChannelList );
      (*theInelastic[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
      G4int itry = 0;
      do
      {
	(*theInelastic[i]).Register(&theNFS, "F01"); // has
	(*theInelastic[i]).Register(&theNXFS, "F02");
	(*theInelastic[i]).Register(&the2NDFS, "F03");
 	(*theInelastic[i]).Register(&the2NFS, "F04"); // has, E Done
 	(*theInelastic[i]).Register(&the3NFS, "F05"); // has, E Done
  	(*theInelastic[i]).Register(&theNAFS, "F06");
	(*theInelastic[i]).Register(&theN3AFS, "F07");
	(*theInelastic[i]).Register(&the2NAFS, "F08");
	(*theInelastic[i]).Register(&the3NAFS, "F09");
	(*theInelastic[i]).Register(&theNPFS, "F10");
	(*theInelastic[i]).Register(&theN2AFS, "F11");
	(*theInelastic[i]).Register(&the2N2AFS, "F12");
	(*theInelastic[i]).Register(&theNDFS, "F13");
	(*theInelastic[i]).Register(&theNTFS, "F14");
	(*theInelastic[i]).Register(&theNHe3FS, "F15");
	(*theInelastic[i]).Register(&theND2AFS, "F16");
	(*theInelastic[i]).Register(&theNT2AFS, "F17");
	(*theInelastic[i]).Register(&the4NFS, "F18"); // has, E Done
	(*theInelastic[i]).Register(&the2NPFS, "F19");
	(*theInelastic[i]).Register(&the3NPFS, "F20");
	(*theInelastic[i]).Register(&theN2PFS, "F21");
	(*theInelastic[i]).Register(&theNPAFS, "F22");
  	(*theInelastic[i]).Register(&thePFS, "F23");
	(*theInelastic[i]).Register(&theDFS, "F24");
	(*theInelastic[i]).Register(&theTFS, "F25");
	(*theInelastic[i]).Register(&theHe3FS, "F26");
	(*theInelastic[i]).Register(&theAFS, "F27");
	(*theInelastic[i]).Register(&the2AFS, "F28");
	(*theInelastic[i]).Register(&the3AFS, "F29");
	(*theInelastic[i]).Register(&the2PFS, "F30");
	(*theInelastic[i]).Register(&thePAFS, "F31");
	(*theInelastic[i]).Register(&theD2AFS, "F32");
	(*theInelastic[i]).Register(&theT2AFS, "F33");
	(*theInelastic[i]).Register(&thePDFS, "F34");
	(*theInelastic[i]).Register(&thePTFS, "F35");
	(*theInelastic[i]).Register(&theDAFS, "F36");
	(*theInelastic[i]).RestartRegistration();
        itry++;
      }
      while( !(*theInelastic[i]).HasDataInAnyFinalState() && itry < 6 );
                                                              // 6 is corresponding to the value(5) of G4NeutronHPChannel. TK  

      if ( itry == 6 ) 
      {
         // No Final State at all.
         G4bool exceptional = false;
         if ( (*(G4Element::GetElementTable()))[i]->GetNumberOfIsotopes() == 1 )
         {
            if ( (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetZ() == 1 && (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetN() == 1 ) exceptional = true;  //1H
         } 
         if ( !exceptional ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this element");
      }

    }
*/
  }

  G4NeutronHPInelastic::~G4NeutronHPInelastic()
  {
//    delete [] theInelastic;
     for ( std::vector<G4NeutronHPChannelList*>::iterator 
           it = (*theInelastic).begin() ; it != (*theInelastic).end() ; it++ )
     {
        delete *it;
     }
     (*theInelastic).clear();
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  
  G4HadFinalState * G4NeutronHPInelastic::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus )
  {
    //if ( numEle < (G4int)G4Element::GetNumberOfElements() ) addChannelForNewElement();
    G4NeutronHPManager::GetInstance()->OpenReactionWhiteBoard();
    const G4Material * theMaterial = aTrack.GetMaterial();
    G4int n = theMaterial->GetNumberOfElements();
    G4int index = theMaterial->GetElement(0)->GetIndex();
    G4int it=0;
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
        //xSec[i] = theInelastic[index].GetXsec(aThermalE.GetThermalEnergy(aTrack,
        xSec[i] = (*(*theInelastic)[index]).GetXsec(aThermalE.GetThermalEnergy(aTrack,
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
        it = i;
        //if(random<=running/sum) break;
        if( sum == 0 || random<=running/sum) break;
      }
      delete [] xSec;
    }

    //return theInelastic[index].ApplyYourself(theMaterial->GetElement(it), aTrack);
    G4HadFinalState* result = (*(*theInelastic)[index]).ApplyYourself(theMaterial->GetElement(it), aTrack);

    //Overwrite target parameters
    aNucleus.SetParameters(G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
    const G4Element* target_element = (*G4Element::GetElementTable())[index];
    const G4Isotope* target_isotope=NULL;
    G4int iele = target_element->GetNumberOfIsotopes();
    for ( G4int j = 0 ; j != iele ; j++ ) { 
       target_isotope=target_element->GetIsotope( j );
       if ( target_isotope->GetN() == G4NeutronHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break; 
    }
    //G4cout << "Target Material of this reaction is " << theMaterial->GetName() << G4endl;
    //G4cout << "Target Element of this reaction is " << target_element->GetName() << G4endl;
    //G4cout << "Target Isotope of this reaction is " << target_isotope->GetName() << G4endl;
    aNucleus.SetIsotope( target_isotope );

    G4NeutronHPManager::GetInstance()->CloseReactionWhiteBoard();
    return result;
  }

const std::pair<G4double, G4double> G4NeutronHPInelastic::GetFatalEnergyCheckLevels() const
{
      // max energy non-conservation is mass of heavy nucleus
//      if ( getenv("G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE") ) return std::pair<G4double, G4double>(5*perCent,250*GeV);
      // This should be same to the hadron default value
//      return std::pair<G4double, G4double>(10*perCent,10*GeV);
      return std::pair<G4double, G4double>(10*perCent,DBL_MAX);
}

/*
void G4NeutronHPInelastic::addChannelForNewElement()
{
   for ( G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++ ) 
   {
      G4cout << "G4NeutronHPInelastic Prepairing Data for the new element of " << (*(G4Element::GetElementTable()))[i]->GetName() << G4endl;

      theInelastic.push_back( new G4NeutronHPChannelList );
      (*theInelastic[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
      G4int itry = 0;
      do
      {
         (*theInelastic[i]).Register(&theNFS, "F01"); // has
         (*theInelastic[i]).Register(&theNXFS, "F02");
         (*theInelastic[i]).Register(&the2NDFS, "F03");
         (*theInelastic[i]).Register(&the2NFS, "F04"); // has, E Done
         (*theInelastic[i]).Register(&the3NFS, "F05"); // has, E Done
         (*theInelastic[i]).Register(&theNAFS, "F06");
         (*theInelastic[i]).Register(&theN3AFS, "F07");
         (*theInelastic[i]).Register(&the2NAFS, "F08");
         (*theInelastic[i]).Register(&the3NAFS, "F09");
         (*theInelastic[i]).Register(&theNPFS, "F10");
         (*theInelastic[i]).Register(&theN2AFS, "F11");
         (*theInelastic[i]).Register(&the2N2AFS, "F12");
         (*theInelastic[i]).Register(&theNDFS, "F13");
         (*theInelastic[i]).Register(&theNTFS, "F14");
         (*theInelastic[i]).Register(&theNHe3FS, "F15");
         (*theInelastic[i]).Register(&theND2AFS, "F16");
         (*theInelastic[i]).Register(&theNT2AFS, "F17");
         (*theInelastic[i]).Register(&the4NFS, "F18"); // has, E Done
         (*theInelastic[i]).Register(&the2NPFS, "F19");
         (*theInelastic[i]).Register(&the3NPFS, "F20");
         (*theInelastic[i]).Register(&theN2PFS, "F21");
         (*theInelastic[i]).Register(&theNPAFS, "F22");
         (*theInelastic[i]).Register(&thePFS, "F23");
         (*theInelastic[i]).Register(&theDFS, "F24");
         (*theInelastic[i]).Register(&theTFS, "F25");
         (*theInelastic[i]).Register(&theHe3FS, "F26");
         (*theInelastic[i]).Register(&theAFS, "F27");
         (*theInelastic[i]).Register(&the2AFS, "F28");
         (*theInelastic[i]).Register(&the3AFS, "F29");
         (*theInelastic[i]).Register(&the2PFS, "F30");
         (*theInelastic[i]).Register(&thePAFS, "F31");
         (*theInelastic[i]).Register(&theD2AFS, "F32");
         (*theInelastic[i]).Register(&theT2AFS, "F33");
         (*theInelastic[i]).Register(&thePDFS, "F34");
         (*theInelastic[i]).Register(&thePTFS, "F35");
         (*theInelastic[i]).Register(&theDAFS, "F36");
         (*theInelastic[i]).RestartRegistration();
         itry++;
      }
      while( !(*theInelastic[i]).HasDataInAnyFinalState() && itry < 6 );
                                                                  // 6 is corresponding to the value(5) of G4NeutronHPChannel. TK  

      if ( itry == 6 ) 
      {
         // No Final State at all.
         G4bool exceptional = false;
         if ( (*(G4Element::GetElementTable()))[i]->GetNumberOfIsotopes() == 1 )
         {
            if ( (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetZ() == 1 && (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetN() == 1 ) exceptional = true;  //1H
         } 
         if ( !exceptional ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this element");
      }
   }

   numEle = (G4int)G4Element::GetNumberOfElements();
}
*/

G4int G4NeutronHPInelastic::GetVerboseLevel() const
{
   return G4NeutronHPManager::GetInstance()->GetVerboseLevel();
}
void G4NeutronHPInelastic::SetVerboseLevel( G4int newValue ) 
{
   G4NeutronHPManager::GetInstance()->SetVerboseLevel(newValue);
}
void G4NeutronHPInelastic::BuildPhysicsTable(const G4ParticleDefinition&)
{
   G4NeutronHPManager* hpmanager = G4NeutronHPManager::GetInstance();

   theInelastic = hpmanager->GetInelasticFinalStates();

   if ( !G4Threading::IsWorkerThread() ) {

      if ( theInelastic == NULL ) theInelastic = new std::vector<G4NeutronHPChannelList*>;

      if ( numEle == (G4int)G4Element::GetNumberOfElements() ) return;

      if ( theInelastic->size() == G4Element::GetNumberOfElements() ) {
         numEle = G4Element::GetNumberOfElements();
         return;
      }

      if ( !getenv("G4NEUTRONHPDATA") ) 
         throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
      dirName = getenv("G4NEUTRONHPDATA");
      G4String tString = "/Inelastic";
      dirName = dirName + tString;

      for (G4int i = numEle ; i < (G4int)G4Element::GetNumberOfElements() ; i++)
      { 
         (*theInelastic).push_back( new G4NeutronHPChannelList );
         (*(*theInelastic)[i]).Init((*(G4Element::GetElementTable()))[i], dirName);
         G4int itry = 0;
         do
         {
	(*(*theInelastic)[i]).Register(&theNFS,    "F01"); // has
	(*(*theInelastic)[i]).Register(&theNXFS,    "F02");
	(*(*theInelastic)[i]).Register(&the2NDFS,    "F03");
    	(*(*theInelastic)[i]).Register(&the2NFS, "F04"); // has, E Done
    	(*(*theInelastic)[i]).Register(&the3NFS, "F05"); // has, E Done
     	(*(*theInelastic)[i]).Register(&theNAFS, "F06");
	(*(*theInelastic)[i]).Register(&theN3AFS,    "F07");
	(*(*theInelastic)[i]).Register(&the2NAFS,    "F08");
	(*(*theInelastic)[i]).Register(&the3NAFS,    "F09");
	(*(*theInelastic)[i]).Register(&theNPFS,    "F10");
	(*(*theInelastic)[i]).Register(&theN2AFS,    "F11");
	(*(*theInelastic)[i]).Register(&the2N2AFS,    "F12");
	(*(*theInelastic)[i]).Register(&theNDFS,    "F13");
	(*(*theInelastic)[i]).Register(&theNTFS,    "F14");
	(*(*theInelastic)[i]).Register(&theNHe3FS,    "F15");
	(*(*theInelastic)[i]).Register(&theND2AFS,    "F16");
	(*(*theInelastic)[i]).Register(&theNT2AFS,    "F17");
	(*(*theInelastic)[i]).Register(&the4NFS,    "F18"); // has, E Done
	(*(*theInelastic)[i]).Register(&the2NPFS,    "F19");
	(*(*theInelastic)[i]).Register(&the3NPFS,    "F20");
	(*(*theInelastic)[i]).Register(&theN2PFS,    "F21");
	(*(*theInelastic)[i]).Register(&theNPAFS,    "F22");
     	(*(*theInelastic)[i]).Register(&thePFS, "F23");
	(*(*theInelastic)[i]).Register(&theDFS,    "F24");
	(*(*theInelastic)[i]).Register(&theTFS,    "F25");
	(*(*theInelastic)[i]).Register(&theHe3FS,    "F26");
	(*(*theInelastic)[i]).Register(&theAFS,    "F27");
	(*(*theInelastic)[i]).Register(&the2AFS,    "F28");
	(*(*theInelastic)[i]).Register(&the3AFS,    "F29");
	(*(*theInelastic)[i]).Register(&the2PFS,    "F30");
	(*(*theInelastic)[i]).Register(&thePAFS,    "F31");
	(*(*theInelastic)[i]).Register(&theD2AFS,    "F32");
	(*(*theInelastic)[i]).Register(&theT2AFS,    "F33");
	(*(*theInelastic)[i]).Register(&thePDFS,    "F34");
	(*(*theInelastic)[i]).Register(&thePTFS,    "F35");
	(*(*theInelastic)[i]).Register(&theDAFS,    "F36");
	(*(*theInelastic)[i]).RestartRegistration();
           itry++;
         }
         while( !(*(*theInelastic)[i]).HasDataInAnyFinalState() && itry < 6 );
                                                                 // 6 is corresponding to the value(5) of G4NeutronHPChannel. TK  

         if ( itry == 6 ) 
         {
            // No Final State at all.
            G4bool exceptional = false;
            if ( (*(G4Element::GetElementTable()))[i]->GetNumberOfIsotopes() == 1 )
            {
               if ( (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetZ() == 1 && (*(G4Element::GetElementTable()))[i]->GetIsotope( 0 )->GetN() == 1 ) exceptional = true;  //1H
            } 
            if ( !exceptional ) throw G4HadronicException(__FILE__, __LINE__, "Channel: Do not know what to do with this element");
         }

      }
      hpmanager->RegisterInelasticFinalStates( theInelastic );

   }
   numEle = G4Element::GetNumberOfElements();
}
