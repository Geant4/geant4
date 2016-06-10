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
// $Id: G4ChipsElasticModel.cc 90228 2015-05-21 08:49:57Z gcosmo $
//
//---------------------------------------------------------------------
//
// Geant4 Class : G4ChipsElasticModel
//
// Author : V.Ivanchenko 29 June 2009 
//  
// Modified:
// 13.01.10: M.Kosov: Use G4Q(Pr/Neut)ElasticCS instead of G4QElasticCS
//
//---------------------------------------------------------------------
// CHIPS model of hadron elastic scattering
//

#include "G4ChipsElasticModel.hh"
#include "G4ParticleDefinition.hh"
#include <iostream>

#include "G4CrossSectionDataSetRegistry.hh"

G4ChipsElasticModel::G4ChipsElasticModel() : G4HadronElastic("hElasticCHIPS")
{
    pxsManager    = (G4ChipsProtonElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name());
    nxsManager    = (G4ChipsNeutronElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name());
    PBARxsManager = (G4ChipsAntiBaryonElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsAntiBaryonElasticXS::Default_Name());
    PIPxsManager  = (G4ChipsPionPlusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsPionPlusElasticXS::Default_Name());
    PIMxsManager  = (G4ChipsPionMinusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsPionMinusElasticXS::Default_Name());
    KPxsManager   = (G4ChipsKaonPlusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusElasticXS::Default_Name());
    KMxsManager   = (G4ChipsKaonMinusElasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusElasticXS::Default_Name());
  //Description();
}

G4ChipsElasticModel::~G4ChipsElasticModel()
{}

void G4ChipsElasticModel::ModelDescription(std::ostream& outFile) const
{

    outFile << "The G4ChipsElasticModel model performs hadron-nucleus elastic\n"
            << "scattering using the parameterized elastic cross sections\n"
            << "of M. Kossov\n";

}


G4double 
G4ChipsElasticModel::SampleInvariantT(const G4ParticleDefinition* p, 
				 G4double plab, G4int Z, G4int A)
{
  G4int N = A - Z;
  if(Z == 1 && N == 2)      { N = 1; }
  else if(Z == 2 && N == 1) { N = 2; }
  G4int projPDG = p->GetPDGEncoding();
  G4double cs = 0.;
  if     (projPDG==2212) { cs = pxsManager->GetChipsCrossSection(plab,Z,N,projPDG); }
  else if(projPDG==2112) { cs = nxsManager->GetChipsCrossSection(plab,Z,N,projPDG); }
  else if(projPDG==-2212){ cs = PBARxsManager->GetChipsCrossSection(plab,Z,N,projPDG); } //Pbar
  else if(projPDG== 211) { cs = PIPxsManager->GetChipsCrossSection(plab,Z,N,projPDG); } // Pi+
  else if(projPDG==-211) { cs = PIMxsManager->GetChipsCrossSection(plab,Z,N,projPDG); } // Pi-
  else if(projPDG== 321) { cs = KPxsManager->GetChipsCrossSection(plab,Z,N,projPDG); } // K+
  else if(projPDG==-321) { cs = KMxsManager->GetChipsCrossSection(plab,Z,N,projPDG); } // K-

  G4double t = 0.0;
  if(cs > 0.0)
  {
    if     (projPDG== 2212) { t = pxsManager->GetExchangeT(Z,N,projPDG); }
    else if(projPDG== 2112) { t = nxsManager->GetExchangeT(Z,N,projPDG); }
    else if(projPDG==-2212) { t = PBARxsManager->GetExchangeT(Z,N,projPDG); } // Pbar
    else if(projPDG==  211) { t = PIPxsManager->GetExchangeT(Z,N,projPDG); }  // Pi+
    else if(projPDG== -211) { t = PIMxsManager->GetExchangeT(Z,N,projPDG); }  // Pi-
    else if(projPDG==  321) { t = KPxsManager->GetExchangeT(Z,N,projPDG); }   // K+
    else if(projPDG== -321) { t = KMxsManager->GetExchangeT(Z,N,projPDG); }   // K-
  }
  else  { t = G4HadronElastic::SampleInvariantT(p, plab, Z, A); }
  return t;
}

