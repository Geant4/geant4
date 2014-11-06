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
//
#ifndef G4NeutronHPIsoData_h
#define G4NeutronHPIsoData_h 1

 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Has the Cross-section data for on isotope.
 
#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
// #include <strstream>
#include <stdlib.h>
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPNames.hh"

class G4NeutronHPIsoData
{
public:

  G4NeutronHPIsoData()
  {
    theChannelData = 0;
    theFissionData = 0;
    theCaptureData = 0;
    theElasticData = 0;
    theInelasticData = 0;
  }
  
  ~G4NeutronHPIsoData(){if(theChannelData!=0) delete theChannelData;}
  
  inline G4double GetXsec(G4double energy)
  {
    return std::max(0., theChannelData->GetXsec(energy));
  }

  //G4bool Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType);
  G4bool Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType){ G4int M = 0 ; return Init( A, Z, M, abun, dirName, aFSType); };
  G4bool Init(G4int A, G4int Z, G4int M, G4double abun, G4String dirName, G4String aFSType);
  
  //void Init(G4int A, G4int Z, G4double abun); //fill PhysicsVector for this Isotope
  void Init(G4int A, G4int Z, G4double abun)  { G4int M =0;
  Init( A, Z, M, abun); }; 
  void Init(G4int A, G4int Z, G4int M, G4double abun); //fill PhysicsVector for this Isotope
  
  inline G4NeutronHPVector * MakeElasticData()
    {return theElasticData;}
  inline G4NeutronHPVector * MakeFissionData()
    {return theFissionData;}
  inline G4NeutronHPVector * MakeCaptureData()
    {return theCaptureData;}
  inline G4NeutronHPVector * MakeInelasticData()
    {return theInelasticData;}
  inline G4NeutronHPVector * MakeChannelData()
    {return theChannelData;}

  G4String GetName(G4int A, G4int Z, G4String base, G4String rest);
  
  inline void FillChannelData(G4NeutronHPVector * aBuffer)
  {
    if(theChannelData!=0) throw G4HadronicException(__FILE__, __LINE__, "IsoData has channel full already!!!");
    theChannelData = new G4NeutronHPVector;
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

  G4NeutronHPVector * theFissionData;
  G4NeutronHPVector * theCaptureData;
  G4NeutronHPVector * theElasticData;
  G4NeutronHPVector * theInelasticData;
  G4NeutronHPVector * theChannelData;

  G4String theFileName;
  G4NeutronHPNames theNames;
};

#endif
