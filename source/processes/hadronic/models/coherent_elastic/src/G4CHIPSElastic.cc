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
// $Id$
//
//---------------------------------------------------------------------
//
// Geant4 Class : G4CHIPSElastic
//
// Author : V.Ivanchenko 29 June 2009 
//  
// Modified:
// 13.01.10: M.Kosov: Use G4Q(Pr/Neut)ElasticCS instead of G4QElasticCS
//
//---------------------------------------------------------------------
// CHIPS model of hadron elastic scattering
//

#include "G4CHIPSElastic.hh"
#include "G4VQCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"

#include "G4QAntiBaryonElasticCrossSection.hh"        // Uzhi
#include "G4QPionPlusElasticCrossSection.hh"          // Uzhi
#include "G4QPionMinusElasticCrossSection.hh"         // Uzhi
#include "G4QKaonPlusElasticCrossSection.hh"          // Uzhi
#include "G4QKaonMinusElasticCrossSection.hh"         // Uzhi
#include <iostream>


G4VQCrossSection* G4CHIPSElastic::pxsManager = 0;
G4VQCrossSection* G4CHIPSElastic::nxsManager = 0;

G4VQCrossSection* G4CHIPSElastic::PBARxsManager = 0;   // Uzhi
G4VQCrossSection* G4CHIPSElastic::PIPxsManager = 0;
G4VQCrossSection* G4CHIPSElastic::PIMxsManager = 0;
G4VQCrossSection* G4CHIPSElastic::KPxsManager = 0;
G4VQCrossSection* G4CHIPSElastic::KMxsManager = 0;

G4CHIPSElastic::G4CHIPSElastic() : G4HadronElastic("hElasticCHIPS")
{
  if(!pxsManager)
  {
    pxsManager    = G4QProtonElasticCrossSection::GetPointer();
    nxsManager    = G4QNeutronElasticCrossSection::GetPointer();

    PBARxsManager = G4QAntiBaryonElasticCrossSection::GetPointer(); // Uzhi
    PIPxsManager  = G4QPionPlusElasticCrossSection::GetPointer();   // Uzhi
    PIMxsManager  = G4QPionMinusElasticCrossSection::GetPointer();  // Uzhi
    KPxsManager   = G4QKaonPlusElasticCrossSection::GetPointer();   // Uzhi
    KMxsManager   = G4QKaonMinusElasticCrossSection::GetPointer();  // Uzhi
  }
  //Description();
}

G4CHIPSElastic::~G4CHIPSElastic()
{}

void G4CHIPSElastic::Description() const
{
  char* dirName = getenv("G4PhysListDocDir");
  if (dirName) {
    std::ofstream outFile;
    G4String outFileName = GetModelName() + ".html";
    G4String pathName = G4String(dirName) + "/" + outFileName;
    outFile.open(pathName);
    outFile << "<html>\n";
    outFile << "<head>\n";

    outFile << "<title>Description of G4CHIPSElastic</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";

    outFile << "The G4CHIPSElastic model performs hadron-nucleus elastic\n"
            << "scattering using the parameterized elastic cross sections\n"
            << "of M. Kossov\n";

    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}


G4double 
G4CHIPSElastic::SampleInvariantT(const G4ParticleDefinition* p, 
				 G4double plab, G4int Z, G4int A)
{
  G4int N = A - Z;
  if(Z == 1 && N == 2)      { N = 1; }
  else if(Z == 2 && N == 1) { N = 2; }
  G4int projPDG = p->GetPDGEncoding();
  G4double cs = 0.;
  if     (projPDG==2212) { cs = pxsManager->GetCrossSection(false,plab,Z,N,projPDG); }
  else if(projPDG==2112) { cs = nxsManager->GetCrossSection(false,plab,Z,N,projPDG); }
  else if(projPDG==-2212){ cs = PBARxsManager->GetCrossSection(false,plab,Z,N,projPDG); } //Pbar
  else if(projPDG== 211) { cs = PIPxsManager->GetCrossSection(false,plab,Z,N,projPDG); } // Pi+
  else if(projPDG==-211) { cs = PIMxsManager->GetCrossSection(false,plab,Z,N,projPDG); } // Pi-
  else if(projPDG== 321) { cs = KPxsManager->GetCrossSection(false,plab,Z,N,projPDG); } // K+
  else if(projPDG==-321) { cs = KMxsManager->GetCrossSection(false,plab,Z,N,projPDG); } // K-

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

