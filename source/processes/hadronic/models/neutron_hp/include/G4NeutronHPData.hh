// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPData.hh,v 1.1 1999-01-07 16:12:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Has the Cross-section data for all materials.
 
#ifndef G4NeutronHPData_h
#define G4NeutronHPData_h 1
#include "globals.hh"
#include "G4Element.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPElementData.hh"

class G4NeutronHPData
{
public:

  G4NeutronHPData();

  ~G4NeutronHPData();

  inline G4PhysicsVector * MakePhysicsVector(G4Element * thE, G4NeutronHPFissionData * theP)
  {
     return DoPhysicsVector(theData[thE->GetIndex()].GetData(theP));
  }
  inline G4PhysicsVector * MakePhysicsVector(G4Element * thE, G4NeutronHPCaptureData * theP)
  {
     return DoPhysicsVector(theData[thE->GetIndex()].GetData(theP));
  }
  inline G4PhysicsVector * MakePhysicsVector(G4Element * thE, G4NeutronHPElasticData * theP)
  {
     return DoPhysicsVector(theData[thE->GetIndex()].GetData(theP));
  }
  inline G4PhysicsVector * MakePhysicsVector(G4Element * thE, G4NeutronHPInelasticData * theP)
  {
//     G4cout << "entered G4NeutronHPData::MakePhysicsVector!!!"<<endl;
//     G4cout << "thE->GetIndex()="<<thE->GetIndex()<<endl;
     return DoPhysicsVector(theData[thE->GetIndex()].GetData(theP));
  }

  G4PhysicsVector * DoPhysicsVector(G4NeutronHPVector * theVector);
  
  static G4NeutronHPData * Instance()
  {
    if(theCrossSectionData==NULL) theCrossSectionData = new G4NeutronHPData;
    return theCrossSectionData;
  }
  
private:

  G4NeutronHPElementData * theData;
  G4int numEle;
  
  static G4NeutronHPData * theCrossSectionData;

};

#endif
