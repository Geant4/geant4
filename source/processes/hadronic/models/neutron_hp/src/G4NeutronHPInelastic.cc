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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
// $Id: G4NeutronHPInelastic.cc,v 1.12 2001-10-19 11:54:51 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4NeutronHPInelastic.hh"

  G4NeutronHPInelastic::G4NeutronHPInelastic()
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 20.*MeV );
    system("echo $NeutronHPCrossSections");
//    G4cout << " entering G4NeutronHPInelastic constructor"<<G4endl;
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    dirName = getenv("NeutronHPCrossSections");
    G4String tString = "/Inelastic/";
    dirName = dirName + tString;
    numEle = G4Element::GetNumberOfElements();
    theInelastic = new G4NeutronHPChannelList[numEle];

    for (G4int i=0; i<numEle; i++)
    { 
      theInelastic[i].Init((*(G4Element::GetElementTable()))[i], dirName);
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
      }
      while(!theInelastic[i].HasDataInAnyFinalState());
    }
  }

  G4NeutronHPInelastic::~G4NeutronHPInelastic()
  {
    delete [] theInelastic;
  }
  
  #include "G4NeutronHPThermalBoost.hh"
  
  G4VParticleChange * G4NeutronHPInelastic::ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
  {
    G4Material * theMaterial = aTrack.GetMaterial();
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
        xSec[i] = theInelastic[index].GetXsec(aThermalE.GetThermalEnergy(aTrack.GetDynamicParticle(),
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
        if(random<=running/sum) break;
      }
      delete [] xSec;
    }
    return theInelastic[index].ApplyYourself(theMaterial->GetElement(it), aTrack);
  }
