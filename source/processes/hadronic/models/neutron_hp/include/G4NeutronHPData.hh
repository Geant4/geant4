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
// $Id: G4NeutronHPData.hh,v 1.4 2001-07-11 10:06:57 gunter Exp $
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
//     G4cout << "entered G4NeutronHPData::MakePhysicsVector!!!"<<G4endl;
//     G4cout << "thE->GetIndex()="<<thE->GetIndex()<<G4endl;
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
