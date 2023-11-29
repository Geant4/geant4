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
/// \file scavenger/src/ParserChemReaction.cc
/// \brief Implementation of the scavenger::ParserChemReaction class

#include "ParserChemReaction.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ParserChemReaction::ReplaceString(G4String &aString, const G4String &from,
                                       const G4String &to) {
  if (G4StrUtil::contains(aString, from)) {
    size_t startPosition = 0;
    while ((startPosition = aString.find(from, startPosition))
           != std::string::npos) {
      aString.replace(startPosition, from.length(), to);
      startPosition += to.length();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParserChemReaction::ImplementReaction(const G4String &reactant1,
                                           const G4String &reactant2,
                                           const std::vector<G4String> &product,
                                           const G4double &reactionRate,
                                           const G4String &type) {
  if (type == "I") {
    fListReactant1[0].push_back(reactant1);
    fListReactant2[0].push_back(reactant2);
    fListProduct[0].push_back(product);
    fListRate[0].push_back(reactionRate);
  } else if (type == "II") {
    fListReactant1[1].push_back(reactant1);
    fListReactant2[1].push_back(reactant2);
    fListProduct[1].push_back(product);
    fListRate[1].push_back(reactionRate);
  } else if (type == "III") {
    fListReactant1[2].push_back(reactant1);
    fListReactant2[2].push_back(reactant2);
    fListProduct[2].push_back(product);
    fListRate[2].push_back(reactionRate);
  } else if (type == "IV") {
    fListReactant1[3].push_back(reactant1);
    fListReactant2[3].push_back(reactant2);
    fListProduct[3].push_back(product);
    fListRate[3].push_back(reactionRate);
  } else if (type == "VI") {
    fListReactant1[4].push_back(reactant1);
    fListReactant2[4].push_back(reactant2);
    fListProduct[4].push_back(product);
    fListRate[4].push_back(reactionRate);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParserChemReaction::ReadReaction(const G4String &reactionString,
                                      std::vector<G4String> &reactant,
                                      std::vector<G4String> &product,
                                      G4double &reactionRate) {
  reactant.clear();
  product.clear();

  G4bool readReaction = true;
  G4bool readReactant = true;
  G4bool readProduct = false;
  G4bool readRate = false;

  std::stringstream aStream;
  aStream << reactionString;

  while (!aStream.eof() && readReaction) {
    G4String aString;
    aStream >> aString;

    if (G4StrUtil::contains(aString,"#")) {
      readReaction = false;
    } else if (readReactant) {
      if (aString == G4String("->")) {
        readReactant = false;
        readProduct = true;
      } else if (aString != G4String("+")) {
        ReplaceString(aString, G4String("+"), G4String("p"));
        ReplaceString(aString, G4String("-"), G4String("m"));

        if (reactant.size() < 2) {
          reactant.push_back(aString);
        }
      }
    } else if (readProduct) {
      if (aString == G4String(",")) {
        readProduct = false;
        readRate = true;
      } else if (aString != G4String("+") && !G4StrUtil::contains(aString,"[")
                 && !G4StrUtil::contains(aString,"]")) {
        ReplaceString(aString, G4String("+"), G4String("p"));
        ReplaceString(aString, G4String("-"), G4String("m"));
        product.push_back(aString);
      }
    } else if (readRate) {
      std::stringstream aStreamTmp;
      aStreamTmp << aString;
      aStreamTmp >> reactionRate;

      if (reactant.size() == 1) {
        // For first-order reactions
        reactionRate *= (1 * 1 / s);
      } else {
        reactionRate *= (1e-3 * m3 / (mole * s));
      }

      readRate = false;
      readReaction = false;
    }
  }

  // For first-order reactions
  if (reactant.size() == 1) {
    reactant.emplace_back("NoneM");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double ParserChemReaction::GetScavengerConcentration(const G4String &name) {
  G4double concentration = -1.;

  std::map<G4String, G4double>::iterator it;
  it = fReservoirConcentrationMap.find(name);

  if (it != fReservoirConcentrationMap.end()) {
    concentration = it->second;
  } else {
    G4ExceptionDescription exception;
    exception << "Scavenger is not defined: "
              << "reaction will not be registered!";
    G4Exception("ParserChemReaction::GetScavengerConcentration", "parchem01",
                JustWarning, exception);
  }

  return concentration;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParserChemReaction::ReadReservoir(const G4String &reservoirString) {
  G4double concentration = 0.;
  G4String name = "";

  G4bool readScavenger = true;
  G4bool readName = false;
  G4bool readConcentration = false;

  std::stringstream aStream;
  aStream << reservoirString;

  while (!aStream.eof() && readScavenger) {
    G4String aString;
    aStream >> aString;

    if (G4StrUtil::contains(aString,"#")) {
      readScavenger = false;
    } else if (aString == G4String("scavenger:")) {
      readName = true;
    } else if (readName) {
      name = G4String(aString);
      ReplaceString(name, G4String("+"), G4String("p"));
      ReplaceString(name, G4String("-"), G4String("m"));

      readName = false;
      readConcentration = true;
    } else if (readConcentration) {
      std::stringstream aStreamTmp;
      aStreamTmp << aString;
      aStreamTmp >> concentration;
      concentration *= (mole / (1e-3 * m3));
      readConcentration = false;
      readScavenger = false;
    }
  }

  if (concentration > 0.) {
    if (fReservoirConcentrationMap.count(name) < 1) {
      fReservoirConcentrationMap[name] = concentration;
    } else {
      G4ExceptionDescription exception;
      exception << "Scavenger already defined previously:\n"
                << "scavenger will not be registered!";
      G4Exception("ParserChemReaction::ReadReservoir", "parchem02",
                  JustWarning, exception);
    }
  } else {
    G4ExceptionDescription exception;
    exception << "Null or negative scavenger concentration:\n"
              << "scavenger will not be registered!";
    G4Exception("ParserChemReaction::ReadReservoir", "parchem03",
                JustWarning, exception);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParserChemReaction::AddReaction(const G4String &reactionString,
                                     const G4String &type) {
  std::vector<G4String> reactant;
  std::vector<G4String> product;
  G4double reactionRate = -1;

  G4bool reservoir = false;

  if (type == "VI") {
    reservoir = true;
  }

  ReadReaction(reactionString, reactant, product, reactionRate);

  if (!reactant.empty() && (reactionRate <= 0)) {
    G4ExceptionDescription exception;
    exception << "Null or negative reaction rate: "
              << "reaction will not be registered!";
    G4Exception("ParserChemReaction::AddReaction", "parchem04",
                JustWarning, exception);
    return;
  }

  G4double concentration;

  if (reservoir && (reactant.size() >= 2)) {
    if (G4StrUtil::contains(reactant[0],"[") && G4StrUtil::contains(reactant[0],"]")) {
      ReplaceString(reactant[0], G4String("["), G4String(""));
      ReplaceString(reactant[0], G4String("]"), G4String(""));

      concentration = GetScavengerConcentration(reactant[0]);

      if (concentration != -1) {
        reactionRate *= concentration;
        reactant[0].append("(B)");
        ImplementReaction(reactant[1], reactant[0], product,
                          reactionRate, type);
      }
    } else if (G4StrUtil::contains(reactant[1],"[") && G4StrUtil::contains(reactant[1],"]")) {
      ReplaceString(reactant[1], G4String("["), G4String(""));
      ReplaceString(reactant[1], G4String("]"), G4String(""));

      concentration = GetScavengerConcentration(reactant[1]);

      if (concentration != -1) {
        reactionRate *= concentration;
        reactant[1].append("(B)");
        ImplementReaction(reactant[0], reactant[1], product,
                          reactionRate, type);
      }
    } else if (reactant[1] == "NoneM") {
      // First-order reaction
      ImplementReaction(reactant[0], reactant[1], product, reactionRate, type);
    } else {
      G4ExceptionDescription exception;
      exception << "Missing or unsuitable square brackets:\n"
                << "reaction will not be registered.\n"
                << "Verify the writing of chemical reactions!";
      G4Exception("ParserChemReaction::AddReaction", "parchem05",
                  JustWarning, exception);
    }
  } else if (reactant.size() >= 2) {
    if (!G4StrUtil::contains(reactant[0],"[") && !G4StrUtil::contains(reactant[0],"]")
        && !G4StrUtil::contains(reactant[1],"[") && !G4StrUtil::contains(reactant[1],"]")
        && (reactant[1] != "NoneM")) {
      ImplementReaction(reactant[0], reactant[1], product, reactionRate, type);
    } else if (reactant[1] == "NoneM") {
      G4ExceptionDescription exception;
      exception << "Unsuitable reaction type: "
                << "reaction will not be registered.\n"
                << "For first-order reaction, use reaction type 6.";
      G4Exception("ParserChemReaction::AddReaction", "parchem06",
                  JustWarning, exception);
    } else {
      G4ExceptionDescription exception;
      exception << "Unsuitable square brackets: "
                << "reaction will not be registered.\n"
                << "Verify the writing of chemical reactions!";
      G4Exception("ParserChemReaction::AddReaction", "parchem07",
                  JustWarning, exception);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParserChemReaction::ReadReactionFile(const G4String &fileName) {
  G4String line;
  std::ifstream myFile(fileName);

  if (myFile.is_open()) {
    while (getline(myFile, line)) {
      if (G4StrUtil::contains(line,"type_1")) {
        AddReaction(line, "I");
      } else if (G4StrUtil::contains(line,"type_2")) {
        AddReaction(line, "II");
      } else if (G4StrUtil::contains(line,"type_3")) {
        AddReaction(line, "III");
      } else if (G4StrUtil::contains(line,"type_4")) {
        AddReaction(line, "IV");
      } else if (G4StrUtil::contains(line,"type_6")) {
        AddReaction(line, "VI");
      } else if (G4StrUtil::contains(line,"scavenger:")) {
        ReadReservoir(line);
      } else if (!G4StrUtil::contains(line,"#") && !line.empty()) {
        G4ExceptionDescription exception;
        exception << "Unknown declaration: "
                  << "reaction or scavenger will not be registered.\n"
                  << "Verify the writing of chemical reactions or scavengers!";
        G4Exception("ParserChemReaction::ReadReactionFile", "parchem08",
                    JustWarning, exception);
      }
    }

    myFile.close();
  } else {
    G4ExceptionDescription exception;
    exception << "Chemical reaction file not found.";
    G4Exception("ParserChemReaction::ReadReactionFile", "parchem09",
                JustWarning, exception);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}