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
#include "G4NeutronHPData.hh"
#include "G4LPhysicsFreeVector.hh"

  G4NeutronHPData::G4NeutronHPData()
  {
     numEle = G4Element::GetNumberOfElements();
     for ( G4int i = 0 ; i < numEle ; i++ ) theData.push_back ( new G4NeutronHPElementData );
//     G4cout << "G4NeutronHPData::G4NeutronHPData(): numEle="<<numEle<<G4endl;
     for (G4int i=0; i<numEle; i++)
     {
       (*theData[i]).Init((*(G4Element::GetElementTable()))[i]);
     }
  }
  
  G4NeutronHPData::~G4NeutronHPData()
  {
     for ( std::vector<G4NeutronHPElementData*>::iterator it = theData.begin() ; it != theData.end() ; it++ ) delete *it;
     theData.clear();
  }
  
  G4NeutronHPData * G4NeutronHPData::Instance()
  {
    static G4NeutronHPData theCrossSectionData;
    return &theCrossSectionData;
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

void G4NeutronHPData::addPhysicsVector()
{
   for ( G4int i = numEle; i < (G4int)G4Element::GetNumberOfElements() ; i++ )
   {
      theData.push_back ( new G4NeutronHPElementData );
      (*theData[i]).Init((*(G4Element::GetElementTable()))[i]);
   }
   numEle = G4Element::GetNumberOfElements();
}
