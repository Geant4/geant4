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
#include "globals.hh"
#include "G4ResonanceNames.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"
#include "G4CrossSectionVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4PionPlus.hh"

G4ResonanceNames::G4ResonanceNames()
{
  // Excited N resonances

 nameNstar.push_back("N(1440)+");
 nameNstar.push_back("N(1440)0");
 nameNstar.push_back("N(1520)+");
 nameNstar.push_back("N(1520)0");
 nameNstar.push_back("N(1535)+");
 nameNstar.push_back("N(1535)0");
 nameNstar.push_back("N(1650)+");
 nameNstar.push_back("N(1650)0");
 nameNstar.push_back("N(1675)+");
 nameNstar.push_back("N(1675)0");
 nameNstar.push_back("N(1680)+");
 nameNstar.push_back("N(1680)0");
 nameNstar.push_back("N(1700)+");
 nameNstar.push_back("N(1700)0");
 nameNstar.push_back("N(1710)+");
 nameNstar.push_back("N(1710)0");
 nameNstar.push_back("N(1720)+");
 nameNstar.push_back("N(1720)0");
 nameNstar.push_back("N(1900)+");
 nameNstar.push_back("N(1900)0");
 nameNstar.push_back("N(1990)+");
 nameNstar.push_back("N(1990)0");
 nameNstar.push_back("N(2090)+");
 nameNstar.push_back("N(2090)0");
 nameNstar.push_back("N(2190)+");
 nameNstar.push_back("N(2190)0");
 nameNstar.push_back("N(2220)+");
 nameNstar.push_back("N(2220)0");
 nameNstar.push_back("N(2250)+");
 nameNstar.push_back("N(2250)0");


  // Delta

  G4String d1232Minus("delta-");
  G4String d1232Zero("delta0");
  G4String d1232Plus("delta+");
  G4String d1232PlusPlus("delta++");
  nameDelta.push_back(d1232Minus);
  nameDelta.push_back(d1232Zero);
  nameDelta.push_back(d1232Plus);
  nameDelta.push_back(d1232PlusPlus);


  // Excited Delta resonances

  nameDeltastar.push_back("delta(1600)+");
  nameDeltastar.push_back("delta(1600)++");
  nameDeltastar.push_back("delta(1600)-");
  nameDeltastar.push_back("delta(1600)0");
  nameDeltastar.push_back("delta(1620)+");
  nameDeltastar.push_back("delta(1620)++");
  nameDeltastar.push_back("delta(1620)-");
  nameDeltastar.push_back("delta(1620)0");
  nameDeltastar.push_back("delta(1700)+");
  nameDeltastar.push_back("delta(1700)++");
  nameDeltastar.push_back("delta(1700)-");
  nameDeltastar.push_back("delta(1700)0");
  nameDeltastar.push_back("delta(1900)+");
  nameDeltastar.push_back("delta(1900)++");
  nameDeltastar.push_back("delta(1900)-");
  nameDeltastar.push_back("delta(1900)0");
  nameDeltastar.push_back("delta(1905)+");
  nameDeltastar.push_back("delta(1905)++");
  nameDeltastar.push_back("delta(1905)-");
  nameDeltastar.push_back("delta(1905)0");
  nameDeltastar.push_back("delta(1910)+");
  nameDeltastar.push_back("delta(1910)++");
  nameDeltastar.push_back("delta(1910)-");
  nameDeltastar.push_back("delta(1910)0");
  nameDeltastar.push_back("delta(1920)+");
  nameDeltastar.push_back("delta(1920)++");
  nameDeltastar.push_back("delta(1920)-");
  nameDeltastar.push_back("delta(1920)0");
  nameDeltastar.push_back("delta(1930)+");
  nameDeltastar.push_back("delta(1930)++");
  nameDeltastar.push_back("delta(1930)-");
  nameDeltastar.push_back("delta(1930)0");
  nameDeltastar.push_back("delta(1950)+");
  nameDeltastar.push_back("delta(1950)++");
  nameDeltastar.push_back("delta(1950)-");
  nameDeltastar.push_back("delta(1950)0");
  

  // Lambda 

  nameLambda.push_back("lambda");
  nameLambda.push_back("lambda(1405)");
  nameLambda.push_back("lambda(1520)");
  nameLambda.push_back("lambda(1600)");
  nameLambda.push_back("lambda(1670)");
  nameLambda.push_back("lambda(1690)");
  nameLambda.push_back("lambda(1800)");
  nameLambda.push_back("lambda(1810)");
  nameLambda.push_back("lambda(1820)");
  nameLambda.push_back("lambda(1830)");
  nameLambda.push_back("lambda(1890)");
  nameLambda.push_back("lambda(2100)");
  nameLambda.push_back("lambda(2110)");


  // Sigma 

  nameSigma.push_back("sigma(1385)+");
  nameSigma.push_back("sigma(1385)-");
  nameSigma.push_back("sigma(1385)0");
  nameSigma.push_back("sigma(1660)+");
  nameSigma.push_back("sigma(1660)-");
  nameSigma.push_back("sigma(1660)0");
  nameSigma.push_back("sigma(1670)+");
  nameSigma.push_back("sigma(1670)-");
  nameSigma.push_back("sigma(1670)0");
  nameSigma.push_back("sigma(1750)+");
  nameSigma.push_back("sigma(1750)-");
  nameSigma.push_back("sigma(1750)0");
  nameSigma.push_back("sigma(1775)+");
  nameSigma.push_back("sigma(1775)-");
  nameSigma.push_back("sigma(1775)0");
  nameSigma.push_back("sigma(1915)+");
  nameSigma.push_back("sigma(1915)-");
  nameSigma.push_back("sigma(1915)0");
  nameSigma.push_back("sigma(1940)+");
  nameSigma.push_back("sigma(1940)-");
  nameSigma.push_back("sigma(1940)0");
  nameSigma.push_back("sigma(2030)+");
  nameSigma.push_back("sigma(2030)-");
  nameSigma.push_back("sigma(2030)0");
  

  // Xi

  nameXi.push_back("xi(1530)-");
  nameXi.push_back("xi(1530)0");
  nameXi.push_back("xi(1690)-");
  nameXi.push_back("xi(1690)0");
  nameXi.push_back("xi(1820)-");
  nameXi.push_back("xi(1820)0");
  nameXi.push_back("xi(1950)-");
  nameXi.push_back("xi(1950)0");
  nameXi.push_back("xi(2030)-");
  nameXi.push_back("xi(2030)0");


  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  size_t i;

  // Fill a map with the lowest resonance for each category
  for (i=0; i<nameNstar.size(); i++)
    {
      lowResMap[nameNstar[i]] = particleTable->FindParticle("N(1440)0");
    }

  for (i=0; i<nameDeltastar.size(); i++)
    {
      lowResMap[nameDeltastar[i]] = particleTable->FindParticle("delta0");
    }

  for (i=0; i<nameDelta.size(); i++)
    {
      lowResMap[nameDelta[i]] = particleTable->FindParticle("delta0");
    }

  for (i=0; i<nameLambda.size(); i++)
    {
      lowResMap[nameLambda[i]] = particleTable->FindParticle("lambda");
    }

  for (i=0; i<nameSigma.size(); i++)
    {
      lowResMap[nameSigma[i]] = particleTable->FindParticle("sigma0");
    }

  shortMap["N(1440)0"] = "N(1440)";
  shortMap["N(1440)+"] = "N(1440)";

  shortMap["N(1520)0"] = "N(1520)";;
  shortMap["N(1520)+"] = "N(1520)";

  shortMap["N(1535)0"] = "N(1535)";
  shortMap["N(1535)+"] = "N(1535)";

  shortMap["N(1650)0"] = "N(1650)";
  shortMap["N(1650)+"] = "N(1650)";

  shortMap["N(1675)0"] = "N(1675)";
  shortMap["N(1675)+"] = "N(1675)";

  shortMap["N(1680)0"] = "N(1680)";
  shortMap["N(1680)+"] = "N(1680)";

  shortMap["N(1700)0"] = "N(1700)";
  shortMap["N(1700)+"] = "N(1700)";

  shortMap["N(1710)0"] = "N(1710)";
  shortMap["N(1710)+"] = "N(1710)";

  shortMap["N(1720)0"] = "N(1720)";
  shortMap["N(1720)+"] = "N(1720)";

  shortMap["N(1900)0"] = "N(1900)";
  shortMap["N(1900)+"] = "N(1900)";

  shortMap["N(1990)0"] = "N(1990)";
  shortMap["N(1990)+"] = "N(1990)";

  shortMap["N(2090)0"] = "N(2090)";
  shortMap["N(2090)+"] = "N(2090)";

  shortMap["N(2190)0"] = "N(2190)";
  shortMap["N(2190)+"] = "N(2190)";

  shortMap["N(2220)0"] = "N(2220)";
  shortMap["N(2220)+"] = "N(2220)";

  shortMap["N(2250)0"] = "N(2250)";
  shortMap["N(2250)+"] = "N(2250)";

 
 // Excited Delta

  shortMap["delta(1600)-"] = "delta(1600)";
  shortMap["delta(1600)0"] = "delta(1600)";
  shortMap["delta(1600)+"] = "delta(1600)";
  shortMap["delta(1600)++"] = "delta(1600)";

  shortMap["delta(1620)-"] = "delta(1620)";
  shortMap["delta(1620)0"] = "delta(1620)";
  shortMap["delta(1620)+"] = "delta(1620)";
  shortMap["delta(1620)++"] = "delta(1620)";

  shortMap["delta(1700)-"] = "delta(1700)";
  shortMap["delta(1700)0"] = "delta(1700)";
  shortMap["delta(1700)+"] = "delta(1700)";
  shortMap["delta(1700)++"] = "delta(1700)";

  shortMap["delta(1900)-"] = "delta(1900)";
  shortMap["delta(1900)0"] = "delta(1900)";
  shortMap["delta(1900)+"] = "delta(1900)";
  shortMap["delta(1900)++"] = "delta(1900)";

  shortMap["delta(1905)-"] = "delta(1905)";
  shortMap["delta(1905)0"] = "delta(1905)";
  shortMap["delta(1905)+"] = "delta(1905)";
  shortMap["delta(1905)++"] = "delta(1905)";

  shortMap["delta(1910)-"] = "delta(1910)";
  shortMap["delta(1910)0"] = "delta(1910)";
  shortMap["delta(1910)+"] = "delta(1910)";
  shortMap["delta(1910)++"] = "delta(1910)";

  shortMap["delta(1920)-"] = "delta(1920)";
  shortMap["delta(1920)0"] = "delta(1920)";
  shortMap["delta(1920)+"] = "delta(1920)";
  shortMap["delta(1920)++"] = "delta(1920)";

  shortMap["delta(1930)-"] = "delta(1930)";
  shortMap["delta(1930)0"] = "delta(1930)";
  shortMap["delta(1930)+"] = "delta(1930)";
  shortMap["delta(1930)++"] = "delta(1930)";

  shortMap["delta(1950)-"] = "delta(1950)";
  shortMap["delta(1950)0"] = "delta(1950)";
  shortMap["delta(1950)+"] = "delta(1950)";
  shortMap["delta(1950)++"] = "delta(1950)";

  // Delta

  shortMap["delta-"] = "delta";
  shortMap["delta0"] = "delta";
  shortMap["delta+"] = "delta";
  shortMap["delta++"] = "delta";

}


G4ResonanceNames::~G4ResonanceNames()
{ }


G4bool G4ResonanceNames::operator==(const G4ResonanceNames &right) const
{
  return(this == (G4ResonanceNames*) &right);
}


G4bool G4ResonanceNames::operator!=(const G4ResonanceNames &right) const
{
  return (this != (G4ResonanceNames*) &right);
}


G4double G4ResonanceNames::MinMass(const G4String& name) 
{
  // Cut, from UrQMD (reference still to be detailed)
  static const G4double coeff = 0.001;

  G4double lowMass = 0.;
  
  G4ParticleDefinition* def = 0;

  if (lowResMap.find(name) != lowResMap.end())
    {
      def = lowResMap[name];
    }
  else
    {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      def = particleTable->FindParticle(name);
    }
  if (def != 0)
    {
      G4ParticleDefinition* pion = G4PionPlus::PionPlusDefinition();
      lowMass = (1. + coeff) * def->GetPDGMass() + pion->GetPDGMass();
    }
  else
    { 
      G4cout << "G4ResonanceNames::MinMass - " << name << " not found" << G4endl;
      throw G4HadronicException(__FILE__, __LINE__,  "G4ResonanceNames::MinMass - resonance name not found");
    }
  return lowMass;
}


const G4String G4ResonanceNames::ShortName(const G4String& name) 
{
  G4String shortName = "";
  if (shortMap.find(name) != shortMap.end())
    {
      shortName = shortMap[name];
    }
  return shortName;
}
