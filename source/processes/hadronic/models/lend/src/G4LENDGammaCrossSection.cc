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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4LENDGammaCrossSection.cc                                        //
//  Date:   30 March 2020                                                     //
//  Author: Dennis H. Wright                                                  //
//                                                                            //
//  Description: cross sections for inelastic scattering of gammas from       //
//               nuclei including gamma-induced fission.  This cross section  //
//               is very similar to G4LENDCombinedCrossSection except that    //
//               does not sample elastic or capture reactions since there are //
//               no such data for gammas in GND.                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4LENDGammaCrossSection.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDFissionCrossSection.hh"
#include "Randomize.hh"


G4LENDGammaCrossSection::G4LENDGammaCrossSection(G4ParticleDefinition* pd)
 :G4LENDCrossSection("LENDGammaCrossSection")
{
  proj = pd;
  inelasticXS = new G4LENDInelasticCrossSection(pd);   
  fissionXS = new G4LENDFissionCrossSection(pd);   
}

void G4LENDGammaCrossSection::BuildPhysicsTable(const G4ParticleDefinition& pd)
{
  inelasticXS->BuildPhysicsTable(pd);
  fissionXS->BuildPhysicsTable(pd);
  create_used_target_map();
}

G4double 
G4LENDGammaCrossSection::GetIsoCrossSection(const G4DynamicParticle* dp, 
                                            G4int iZ, G4int iA,
                                            const G4Isotope* isotope,
                                            const G4Element*,
                                            const G4Material* material)
{
  G4double XS = 0.0;
  XS += inelasticXS->GetIsoCrossSection(dp, iZ, iA, isotope, NULL, material);
  XS += fissionXS->GetIsoCrossSection(dp, iZ, iA, isotope, NULL, material);
  //G4cout << "G4LENDGammaCrossSection::GetIsoCrossSection " 
  //       << XS/CLHEP::barn << " [barn]" << G4endl;
  return XS;
}

G4int G4LENDGammaCrossSection::SelectChannel(const G4DynamicParticle* dp,
                                             G4int iZ, G4int iA,
                                             const G4Isotope* isotope,
                                             const G4Element*,
                                             const G4Material* material)
{
  G4int ichannel = -1;
  G4double XSs[2];
  XSs[0] = inelasticXS->GetIsoCrossSection(dp, iZ, iA, isotope, nullptr, material);
  XSs[1] = XSs[0] + fissionXS->GetIsoCrossSection(dp, iZ, iA, isotope, nullptr, material);

  G4double total = XSs[1];

  G4double random = G4UniformRand();
  for (G4int i = 0; i < 2; i++) {
    if (random*total <= XSs[i]) { 
      ichannel = i;
      break; 
    }
  }
   
  return ichannel;
}

