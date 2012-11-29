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
// $Id$
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
//
// 080520 Delete unnecessary dependencies by T. Koi
 
#ifndef G4NeutronHPElementData_h
#define G4NeutronHPElementData_h 1
#include "globals.hh"
#include "G4NeutronHPIsoData.hh"
#include "G4NeutronHPVector.hh"
#include "G4Material.hh"
#include "G4HadronCrossSections.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
//#include "G4NeutronInelasticProcess.hh"
//#include "G4HadronFissionProcess.hh"
//#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4StableIsotopes.hh"

class G4NeutronHPElementData : public G4HadronCrossSections
{
public:

  G4NeutronHPElementData();
  
  ~G4NeutronHPElementData();
  
  void Init(G4Element * theElement);
  
  //void UpdateData(G4int A, G4int Z, G4int index, G4double abundance);
  void UpdateData(G4int A, G4int Z, G4int index, G4double abundance) { G4int M=0; UpdateData( A, Z, M, index, abundance); };
  void UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance);
  
  void Harmonise(G4NeutronHPVector *& theStore, G4NeutronHPVector * theNew);
  
  inline G4NeutronHPVector * GetData(G4NeutronHPFissionData * )
    {return theFissionData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPCaptureData * )
    {return theCaptureData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPElasticData * )
    {return theElasticData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPInelasticData * )
    {return theInelasticData;}

  G4NeutronHPVector * MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPFissionData* theSet);

  G4NeutronHPVector * MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPCaptureData * theSet);

  G4NeutronHPVector * MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPElasticData * theSet);

  G4NeutronHPVector * MakePhysicsVector(G4Element * theElement,
                                      G4ParticleDefinition * theP,
                                      G4NeutronHPInelasticData * theSet);

private:

  G4NeutronHPVector * theFissionData;
  G4NeutronHPVector * theCaptureData;
  G4NeutronHPVector * theElasticData;
  G4NeutronHPVector * theInelasticData;
  G4double precision;
  
  G4NeutronHPVector * theBuffer;
  
  G4NeutronHPIsoData * theIsotopeWiseData;

  G4StableIsotopes theStableOnes;
  
  G4String filename;
  
};

#endif
