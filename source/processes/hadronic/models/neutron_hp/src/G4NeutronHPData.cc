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
//     G4cout << "G4NeutronHPData::G4NeutronHPData(): numEle="<<numEle<<endl;
     for (G4int i=0; i<numEle; i++)
     {
       theData[i].Init((*(G4Element::GetElementTable()))(i));
     }
  }
  
  G4NeutronHPData::~G4NeutronHPData()
  {
  delete [] theData;
  }
  
  G4PhysicsVector * G4NeutronHPData::DoPhysicsVector(G4NeutronHPVector * theVector)
  {
//    G4cout << "Entered G4NeutronHPData::DoPhysicsVector."<<endl;
    G4int len = theVector->GetVectorLength();
//    G4cout <<"zahl der energie-punkte "<< len<<endl;
    if(len==0) return new G4LPhysicsFreeVector(0, 0, 0);
    G4double emin = theVector->GetX(0);
    G4double emax = theVector->GetX(len-1);
//    G4cout <<"zahl der energie-punkte "<< len<<" "<<emin<<" "<<emax<<endl;
    
   //  G4int dummy; cin >> dummy; 
    G4bool flag;
    
    G4LPhysicsFreeVector * theResult = new G4LPhysicsFreeVector(len, emin, emax);
    for (G4int i=0; i<len; i++)
    {
      theResult->PutValues(i, theVector->GetX(i), theVector->GetY(i));
    }
    return theResult;
  }
