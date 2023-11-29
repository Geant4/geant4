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
/*
 * G4ReactionTableMessenger.cc
 *
 *  Created on: Sep 14, 2015
 *      Author: mkaramit
 */

#include <G4ReactionTableMessenger.hh>
#include <G4DNAMolecularReactionTable.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithABool.hh>

//------------------------------------------------------------------------------

G4ReactionTableMessenger::G4ReactionTableMessenger(G4DNAMolecularReactionTable* table)
  : G4UImessenger()
  , fpTable(table)
  , fpActivateReactionUI(new G4UIcmdWithoutParameter("/chem/reaction/UI", this))
{
  fpNewDiffContReaction = new G4UIcmdWithAString("/chem/reaction/new", this);
  fpAddReaction = new G4UIcmdWithAString("/chem/reaction/add", this);
  fpPrintTable = new G4UIcmdWithoutParameter("/chem/reaction/print", this);
}

//------------------------------------------------------------------------------

G4ReactionTableMessenger::~G4ReactionTableMessenger()
{
  if(fpNewDiffContReaction) delete fpNewDiffContReaction;
  if(fpAddReaction) delete fpAddReaction;
  if(fpPrintTable) delete fpPrintTable;
}

//------------------------------------------------------------------------------
void G4ReactionTableMessenger::SetNewValue(G4UIcommand* command,
                                           G4String newValue)
{
  if(command == fpActivateReactionUI.get())
  {
    //assert(false);
      fpTable->Reset();//release reaction data
  }

  if(command == fpNewDiffContReaction)
  {
    std::istringstream iss(newValue);

    G4String species1;
    iss >> species1;

    G4String species2;
    iss >> species2;

    double reactionRate;
    iss >> reactionRate;

//    G4String reactionRateUnit;
//    iss >> reactionRateUnit;

    double dimensionedReactionRate = reactionRate * (1e-3 * m3 / (mole * s));
//        G4UIcmdWithADoubleAndUnit::ConvertToDimensionedDouble((reactionRate
//            + G4String(" ") + reactionRateUnit).c_str());

    auto reactionData =
        new G4DNAMolecularReactionData(dimensionedReactionRate,
                                       species1,
                                       species2);

//    G4String productionRate;
//    iss >> productionRate;
//
//    if(productionRate != "" && productionRate != "X")
//    {
//      double prodRateReal = G4UIcommand::ConvertToDouble(productionRate);
//
//      if(prodRateReal == 0 || isnan(prodRateReal))
//      {
//        G4Exception("G4ReactionTableMessenger",
//                    "WRONG_PRODUCTION_RATE",
//                    FatalException, "");
//      }
//
//      double dimensionedProductionRate = prodRateReal
//          * (1e-3 * m3 / (mole * s));
//      reactionData->SetProductionRate(dimensionedProductionRate);
//    }

    while(iss.eof() == false)
    {
      G4String product;
      iss >> product;

      if(product != "")
      {
        reactionData->AddProduct(product);
      }
      else
      {
        break;
      }
    };

    fpTable->SetReaction(reactionData);
  }
//  else if(command == fpNewPartDiffContReactionByRadius)
//  {
//    std::istringstream iss(newValue);
//
//    G4String species1;
//    iss >> species1;
//
//    G4String species2;
//    iss >> species2;
//
//    double reactionRate;
//    iss >> reactionRate;
//
////    G4String reactionRateUnit;
////    iss >> reactionRateUnit;
//
////    G4String reactionRateUnit;
////    iss >> reactionRateUnit;
//
//    double reactionRadius;
//    iss >> reactionRadius;
//
////    G4String reactionRadiusUnit;
////    iss >> reactionRadiusUnit;
//
//    double dimensionedReactionRate = reactionRate * (1e-3 * m3 / (mole * s));
//
////    double dimensionedReactionRate =
////        G4UIcmdWithADoubleAndUnit::ConvertToDimensionedDouble((reactionRate
////            + " " + reactionRateUnit).c_str());
//
//    double dimensionedReactionRadius = reactionRadius * nm;
////        G4UIcmdWithADoubleAndUnit::ConvertToDimensionedDouble((reactionRadius
////            + " " + reactionRadiusUnit).c_str());
//
//    G4DNAMolecularReactionData* reactionData =
//        new G4DNAMolecularReactionData(dimensionedReactionRate,
//                                       species1,
//                                       species2);
//    reactionData->SetPartiallyDiffusionControlledReaction(dimensionedReactionRate,
//                                                          dimensionedReactionRadius);
//
//    while(iss.eof() == false)
//    {
//      G4String product;
//      iss >> product;
//
//      if(product != "")
//      {
//        reactionData->AddProduct(product);
//      }
//      else
//      {
//        break;
//      }
//    }
//
//    fpTable->SetReaction(reactionData);
//  }
//  else if(command == fpNewPartDiffContReactionByReactionRate)
//  {
//    std::istringstream iss(newValue);
//
//    G4String species1;
//    iss >> species1;
//
//    G4String species2;
//    iss >> species2;
//
//    double reactionRate;
//    iss >> reactionRate;
//
//    //    G4String reactionRateUnit;
//    //    iss >> reactionRateUnit;
//
//    //    G4String reactionRateUnit;
//    //    iss >> reactionRateUnit;
//
//    double activationRate;
//    iss >> activationRate;
//
//    //    G4String reactionRadiusUnit;
//    //    iss >> reactionRadiusUnit;
//
//    double dimensionedReactionRate = reactionRate * (1e-3 * m3 / (mole * s));
//
//    double dimensionedActivationRate = activationRate
//        * (1e-3 * m3 / (mole * s));
//
//    G4DNAMolecularReactionData* reactionData =
//        new G4DNAMolecularReactionData(dimensionedReactionRate,
//                                       species1,
//                                       species2);
//    reactionData->
//      SetPartiallyDiffusionControlledReactionByActivation(dimensionedReactionRate,
//                                                          dimensionedActivationRate);
//
//    while(iss.eof() == false)
//    {
//      G4String product;
//      iss >> product;
//
//      if(product != "")
//      {
//        reactionData->AddProduct(product);
//      }
//      else
//      {
//        break;
//      }
//    }
//
//    fpTable->SetReaction(reactionData);
//  }
  else if(command == fpPrintTable)
  {
    fpTable->PrintTable();
  }
  else if(command == fpAddReaction)
  {
    std::istringstream iss(newValue);

    //--------------------------------------------------------------------------
    // Reactants definition

    G4String species1;
    iss >> species1;

    G4String marker;
    iss >> marker; // peut etre +, ->, |

    G4String species2;

    if(marker == "+")
    {
        iss >> species2;
        iss >> marker; // peut etre ->, |
    }

    //--------------------------------------------------------------------------

    G4DNAMolecularReactionData* reactionData =
                new G4DNAMolecularReactionData(0,
                                               species1,
                                               species2);
    //fpTable->SetReaction(reactionData);

    //--------------------------------------------------------------------------
    // Add products
    if(marker == "->")
    {
      iss >> marker; // doit etre = species name

      while(marker!="|"
          //&& marker!=""
          && iss.eof() == false
          )
      {
        //G4cout << marker << G4endl;
        if(marker == "+")
        {
          iss >> marker; // doit etre species name
          continue;
        }
        if(marker != "H2O")
        {
          reactionData->AddProduct(marker);
        }

        iss >> marker; // peut etre species name, +, |
      };
    }

//    G4cout << "out of 1st loop" << G4endl;

    //--------------------------------------------------------------------------
    // Add reaction rate method
    G4String rateconst_method;
    iss >> rateconst_method;
    if(rateconst_method == "Fix")
    {
      iss >> marker; // must be |
      double reactionRate;
      iss >> reactionRate;

      double dimensionedReactionRate = reactionRate * (1e-3 * m3 / (mole * s));
      reactionData->SetObservedReactionRateConstant(dimensionedReactionRate);
      reactionData->ComputeEffectiveRadius();
      G4String markerType;
      iss >> markerType; // must be |
      if(markerType == "|")
      {
        G4int reactionType;
        iss >> reactionType;
        if(reactionType == 1)
        {
          reactionData->SetReactionType(reactionType);
        }
      }


//      G4String productionRate;
//      iss >> productionRate;
//
//      if(productionRate != "" && productionRate != "X")
//      {
//        double prodRateReal = G4UIcommand::ConvertToDouble(productionRate);
//
//        if(prodRateReal == 0 || isnan(prodRateReal))
//        {
//          G4Exception("G4ReactionTableMessenger",
//                      "WRONG_PRODUCTION_RATE",
//                      FatalException,
//                      "");
//        }
//
//        double dimensionedProductionRate = prodRateReal
//            * (1e-3 * m3 / (mole * s));
//        reactionData->SetProductionRate(dimensionedProductionRate);
//      }
    }
    else if(rateconst_method == "Arr")
    {
      iss >> marker; // must be |
      double A0  = 0;
      double E_R = 0;

      iss >> A0;
      iss >> E_R;
      reactionData->SetArrehniusParameterization(A0, E_R);
    }
    else if(rateconst_method == "Pol")
    {
      iss >> marker; // must be |
      std::vector<double> P = {0, 0, 0, 0, 0};

      size_t i = 0;
      while(i < 4) // could be changed to 5 only if marker is used as delimiter
      {
        double tmp;
        iss >> tmp;
        P[i] = tmp;
//        G4cout << newValue << G4endl;
//        G4cout << tmp << G4endl;
//        G4cout << P[i] << G4endl;
        ++i;
      };
      reactionData->SetPolynomialParameterization(P);
    }
    else if(rateconst_method == "Scale")
    {
      iss >> marker; // must be |
      double temp_K;
      iss >> temp_K;
      double reactionRateCste;
      iss >> reactionRateCste;
      double dimensionedReactionRate = reactionRateCste * (1e-3 * m3 / (mole * s));
      reactionData->SetObservedReactionRateConstant(dimensionedReactionRate);
      reactionData->SetScaledParameterization(temp_K, dimensionedReactionRate);
    }

//    if(iss.eof() == false)
//    {
//      iss >> marker;
//
//      if(marker == "|")
//      {
//        G4String productionRate ;
//        iss >> productionRate;
//
////        G4cout << productionRate << G4endl;
//
//        double dimProductionRate = G4UIcommand::ConvertToDouble(productionRate)* (1e-3 * m3 / (mole * s));
//
//        G4cout << " DIM PROD RATE = " << reactionData->GetReactant1()->GetName()
//             << " + " << reactionData->GetReactant2()->GetName() << " = " << dimProductionRate << G4endl;
//
//        reactionData->SetProductionRate(dimProductionRate);
//      }
//    }
    fpTable->SetReaction(reactionData);
//    G4cout << "Reaction " << species1 << " + " << species2 << " added" << G4endl;
  }
}
