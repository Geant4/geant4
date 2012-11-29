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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4CHIPSElasticXS
//
// Author  Ivantchenko, Geant4, 3-Aug-09
//
// Modifications:
// 31-05-2011 V.Uzhinsky added anti-baryons, Pi+, Pi-, K+, K- cross sections
// 23-08-2011 V.Ivanchenko migration to new design and cleanup
//

#include "G4CHIPSElasticXS.hh"
#include "G4HadronicException.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4VQCrossSection.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"

#include "G4QAntiBaryonElasticCrossSection.hh" // Uzhi
#include "G4QPionMinusElasticCrossSection.hh"  // Uzhi
#include "G4QPionPlusElasticCrossSection.hh"   // Uzhi
#include "G4QKaonMinusElasticCrossSection.hh"  // Uzhi
#include "G4QKaonPlusElasticCrossSection.hh"   // Uzhi

G4CHIPSElasticXS::G4CHIPSElasticXS() 
  :  G4VCrossSectionDataSet("CHIPSElasticXS"),
     theProton(G4Proton::Proton()), 
     theNeutron(G4Neutron::Neutron()),
     thEnergy(19*CLHEP::MeV),
     isInitialized(false)
{
  //  verboseLevel = 0;
  pCManager   = G4QProtonElasticCrossSection::GetPointer();
  nCManager   = G4QNeutronElasticCrossSection::GetPointer();

  PBARxsManager = G4QAntiBaryonElasticCrossSection::GetPointer(); // Uzhi
  PIPxsManager  = G4QPionPlusElasticCrossSection::GetPointer();   // Uzhi
  PIMxsManager  = G4QPionMinusElasticCrossSection::GetPointer();  // Uzhi
  KPxsManager   = G4QKaonPlusElasticCrossSection::GetPointer();   // Uzhi
  KMxsManager   = G4QKaonMinusElasticCrossSection::GetPointer();  // Uzhi
  //Description();
  theParticle   = 0;
}

G4CHIPSElasticXS::~G4CHIPSElasticXS()
{}


void G4CHIPSElasticXS::Description() const
{
  char* dirName = getenv("G4PhysListDocDir");
  if (dirName) {
    std::ofstream outFile;
    G4String outFileName = GetName() + ".html";
    G4String pathName = G4String(dirName) + "/" + outFileName;

    outFile.open(pathName);
    outFile << "<html>\n";
    outFile << "<head>\n";

    outFile << "<title>Description of CHIPS Elastic Cross Section</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";

    outFile << "G4CHIPSElasticXS provides hadron-nuclear elastic scattering\n"
            << "cross sections for protons and neutrons with incident energies\n"
            << "between 19 MeV and X GeV.  These cross sections represent\n"
            << "parameterizations developed by M. Kossov. (more detail)\n";

    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}

G4bool 
G4CHIPSElasticXS::IsIsoApplicable(const G4DynamicParticle* dyn, 
				  G4int Z, G4int /*A*/,
				  const G4Element*, const G4Material*)
{
  return (Z <= 2 && dyn->GetKineticEnergy() > thEnergy);
}

G4double 
G4CHIPSElasticXS::GetIsoCrossSection(const G4DynamicParticle* dyn, 
				     G4int Z, G4int A,
				     const G4Isotope*, const G4Element*, 
				     const G4Material*)
{
  G4int N = A - Z;
  if(Z == 1) {
    if(N > 1) { N = 1; }
  } else if(Z == 2) { N = 2; }

  G4double momentum = dyn->GetTotalMomentum();
  G4int    uPDGcode = dyn->GetPDGcode();
  G4VQCrossSection* CHIPSmanager = 0; 
  G4double cross = 0.0;

  switch(uPDGcode) {
  case 2212:
    CHIPSmanager=pCManager;
    break;
  case 2112:
    CHIPSmanager=nCManager;
    break;
  case -2212:
    CHIPSmanager=PBARxsManager;
    break;
  case -2112:
    CHIPSmanager=PBARxsManager;
    break;
  case 211:
    CHIPSmanager=PIPxsManager;
    break;
  case -211:
    CHIPSmanager=PIMxsManager;
    break;
  case 321:
    CHIPSmanager=KPxsManager;
    break;
  case -321:
    CHIPSmanager=KMxsManager;
    break;
  case 130:
    break;
  case 310:
    break;
  case 311:
    break;
  case -311:
    break;
  default:
    throw G4HadronicException(__FILE__, __LINE__,
			      "G4CHIPSElasticXS: not applicable for a particle"); 
    return cross; 
  }
  if(CHIPSmanager) {
    cross = CHIPSmanager->GetCrossSection(false,momentum,Z,N,uPDGcode);
  } else {
    cross = 0.5*(KPxsManager->GetCrossSection(false,momentum,Z,N,uPDGcode) +
		 KMxsManager->GetCrossSection(false,momentum,Z,N,uPDGcode));
  }
  return cross; 
}
