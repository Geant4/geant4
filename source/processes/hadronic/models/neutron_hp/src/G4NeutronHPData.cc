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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPData.hh"
#include "G4LPhysicsFreeVector.hh"

  G4NeutronHPData * G4NeutronHPData::theCrossSectionData = NULL;

  G4NeutronHPData::G4NeutronHPData()
  {
     numEle = G4Element::GetNumberOfElements();
     theData = new G4NeutronHPElementData[numEle];
//     G4cout << "G4NeutronHPData::G4NeutronHPData(): numEle="<<numEle<<G4endl;
     for (G4int i=0; i<numEle; i++)
     {
       theData[i].Init((*(G4Element::GetElementTable()))[i]);
     }
  }
  
  G4NeutronHPData::~G4NeutronHPData()
  {
  delete [] theData;
  }
  
  G4PhysicsVector * G4NeutronHPData::DoPhysicsVector(G4NeutronHPVector * theVector)
  {
//    G4cout << "Entered G4NeutronHPData::DoPhysicsVector."<<G4endl;
    G4int len = theVector->GetVectorLength();
//    G4cout <<"zahl der energie-punkte "<< len<<G4endl;
    if(len==0) return new G4LPhysicsFreeVector(0, 0, 0);
    G4double emin = theVector->GetX(0);
    G4double emax = theVector->GetX(len-1);
//    G4cout <<"zahl der energie-punkte "<< len<<" "<<emin<<" "<<emax<<G4endl;
    
   //  G4int dummy; G4cin >> dummy;     
    G4LPhysicsFreeVector * theResult = new G4LPhysicsFreeVector(len, emin, emax);
    for (G4int i=0; i<len; i++)
    {
      theResult->PutValues(i, theVector->GetX(i), theVector->GetY(i));
    }
    return theResult;
  }
