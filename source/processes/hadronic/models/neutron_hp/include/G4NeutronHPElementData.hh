// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPElementData.hh,v 1.3 1999-10-22 08:24:57 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPElementData_h
#define G4NeutronHPElementData_h 1
#include "globals.hh"
#include "G4NeutronHPIsoData.hh"
#include "G4NeutronHPVector.hh"
#include "G4Material.hh"
#include "G4HadronCrossSections.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
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
  
  void UpdateData(G4int A, G4int Z, G4int index, G4double abundance);
  
  void Harmonise(G4NeutronHPVector *& theStore, G4NeutronHPVector * theNew);
  
  inline G4NeutronHPVector * GetData(G4NeutronHPFissionData * theP)
    {return theFissionData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPCaptureData * theP)
    {return theCaptureData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPElasticData * theP)
    {return theElasticData;}
  inline G4NeutronHPVector * GetData(G4NeutronHPInelasticData * theP)
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
