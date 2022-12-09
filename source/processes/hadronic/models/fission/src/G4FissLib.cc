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
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
//
//
// neutron_hp -- source file
// J.M. Verbeke, Jan-2007
// A prototype of the low energy neutron transport model.
//
#include "G4FissLib.hh"
#include "G4SystemOfUnits.hh"

G4FissLib::G4FissLib()
 :xSec(0)
{
  SetMinEnergy(0.0);
  SetMaxEnergy(20.*MeV);
  if(!G4FindDataDir("G4NEUTRONHPDATA")) {
     G4cout << "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files." << G4endl;
     throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files.");
  }
  dirName = G4FindDataDir("G4NEUTRONHPDATA");
  G4String tString = "/Fission/";
  dirName = dirName + tString;
  numEle = (G4int)G4Element::GetNumberOfElements();
  theFission = new G4ParticleHPChannel[numEle];

  for (G4int i=0; i<numEle; ++i)
  { 
//    G4cout << "G4FissLib::G4FissLib(): element "<< i << " : " << (*(G4Element::GetElementTable()))[i]->GetZ()<< G4endl;
    if((*(G4Element::GetElementTable()))[i]->GetZ()>89)
    {
      theFission[i].Init((*(G4Element::GetElementTable()))[i], dirName);
      theFission[i].Register(&theLibrary);
    }
  }
}
  
G4FissLib::~G4FissLib()
{
  delete [] theFission;
}
  
G4HadFinalState*
G4FissLib::ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aNucleus)
{
  G4ParticleHPManager::GetInstance()->OpenReactionWhiteBoard();

  const G4Material* theMaterial = aTrack.GetMaterial();
  std::size_t n = theMaterial->GetNumberOfElements();
  std::size_t index = theMaterial->GetElement(0)->GetIndex();

  if (n != 1) {
    xSec = new G4double[n];
    G4double sum = 0;
    G4int i;
    std::size_t imat;
    const G4double * NumAtomsPerVolume = theMaterial->GetVecNbOfAtomsPerVolume();
    G4double rWeight;    
    G4ParticleHPThermalBoost aThermalE;
    for (i = 0; i < (G4int)n; ++i) {
      imat = theMaterial->GetElement(i)->GetIndex();
      rWeight = NumAtomsPerVolume[i];
      xSec[i] = theFission[imat].GetXsec(aThermalE.GetThermalEnergy(aTrack,
		                                                    theMaterial->GetElement(i),
  						      theMaterial->GetTemperature()));
      xSec[i] *= rWeight;
      sum+=xSec[i];
    }

    G4double random = G4UniformRand();
    G4double running = 0;
    for (i = 0; i < (G4int)n; ++i) {
      running += xSec[i];
      index = theMaterial->GetElement(i)->GetIndex();
      if(random<=running/sum) break;
    }
    delete [] xSec;
  }

  //return theFission[index].ApplyYourself(aTrack);
   //Overwrite target parameters
   G4HadFinalState* result = theFission[index].ApplyYourself(aTrack);
   aNucleus.SetParameters(G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA(),G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargZ());
   const G4Element* target_element = (*G4Element::GetElementTable())[index];
   const G4Isotope* target_isotope=NULL;
   G4int iele = (G4int)target_element->GetNumberOfIsotopes();
   for ( G4int j = 0 ; j != iele ; ++j ) {
      target_isotope=target_element->GetIsotope( j );
      if ( target_isotope->GetN() == G4ParticleHPManager::GetInstance()->GetReactionWhiteBoard()->GetTargA() ) break;
   }
   aNucleus.SetIsotope( target_isotope );
   G4ParticleHPManager::GetInstance()->CloseReactionWhiteBoard();
   return result;
}

const std::pair<G4double, G4double> G4FissLib::GetFatalEnergyCheckLevels() const
{
  // max energy non-conservation is mass of heavy nucleus
  return std::pair<G4double, G4double>(10.0*perCent, 350.0*CLHEP::GeV);
}
