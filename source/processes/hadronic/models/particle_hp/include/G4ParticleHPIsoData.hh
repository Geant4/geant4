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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPIsoData_h
#define G4ParticleHPIsoData_h 1

 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Has the Cross-section data for on isotope.
 
#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
// #include <strstream>
#include <stdlib.h>
#include "G4ParticleHPVector.hh"
#include "G4ParticleHPNames.hh"
class G4ParticleDefinition;

class G4ParticleHPIsoData
{
public:

  G4ParticleHPIsoData()
  {
    theChannelData = 0;
    theFissionData = 0;
    theCaptureData = 0;
    theElasticData = 0;
    theInelasticData = 0;
  }
  
  ~G4ParticleHPIsoData(){if(theChannelData!=0) delete theChannelData;}
  
  inline G4double GetXsec(G4double energy)
  {
    return std::max(0., theChannelData->GetXsec(energy));
  }

  //G4bool Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType);
  G4bool Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType){ G4int M = 0 ; return Init( A, Z, M, abun, dirName, aFSType); };
  G4bool Init(G4int A, G4int Z, G4int M, G4double abun, G4String dirName, G4String aFSType);
  
  //void Init(G4int A, G4int Z, G4double abun); //fill PhysicsVector for this Isotope
  void Init(G4int A, G4int Z, G4double abun, G4ParticleDefinition* projectile, const char* dataDirVariable)  { G4int M =0;
    Init( A, Z, M, abun, projectile, dataDirVariable ); }; 
  void Init(G4int A, G4int Z, G4int M, G4double abun, G4ParticleDefinition* projectile, const char* dataDirVariable); //fill PhysicsVector for this Isotope
  
  inline G4ParticleHPVector * MakeElasticData()
    {return theElasticData;}
  inline G4ParticleHPVector * MakeFissionData()
    {return theFissionData;}
  inline G4ParticleHPVector * MakeCaptureData()
    {return theCaptureData;}
  inline G4ParticleHPVector * MakeInelasticData()
    {return theInelasticData;}
  inline G4ParticleHPVector * MakeChannelData()
    {return theChannelData;}

  G4String GetName(G4int A, G4int Z, G4String base, G4String rest);
  
  inline void FillChannelData(G4ParticleHPVector * aBuffer)
  {
    if(theChannelData!=0) throw G4HadronicException(__FILE__, __LINE__, "IsoData has channel full already!!!");
    theChannelData = new G4ParticleHPVector;
    for(G4int i=0; i<aBuffer->GetVectorLength(); i++)
    {
      theChannelData->SetPoint(i, aBuffer->GetPoint(i));
    }
    theChannelData->Hash();
  }
  
  inline void ThinOut(G4double precision)
  {
    if(theFissionData) theFissionData->ThinOut(precision);
    if(theCaptureData) theCaptureData->ThinOut(precision);
    if(theElasticData) theElasticData->ThinOut(precision);
    if(theInelasticData) theInelasticData->ThinOut(precision);
  }
  
private:

  G4ParticleHPVector * theFissionData;
  G4ParticleHPVector * theCaptureData;
  G4ParticleHPVector * theElasticData;
  G4ParticleHPVector * theInelasticData;
  G4ParticleHPVector * theChannelData;

  G4String theFileName;
  G4ParticleHPNames theNames;
};

#endif
