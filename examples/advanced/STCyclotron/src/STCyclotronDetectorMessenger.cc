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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronDetectorMesseger.cc
//

#include "STCyclotronDetectorMessenger.hh"
#include "STCyclotronDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"

/////////////////////////////////////////////////////////////////////////////
STCyclotronDetectorMessenger::STCyclotronDetectorMessenger(STCyclotronDetectorConstruction* detector)
  :fDet(detector)
{
  /////////////////////////////  
  // Change Target parameters//
  /////////////////////////////
    fChangeTarget = new G4UIdirectory("/changeTarget/");
    fChangeTarget -> SetGuidance("Change the Target diameter/thickness/material");

    // change Target diameter
    fChangeTargetDiameterCmd = new G4UIcmdWithADoubleAndUnit("/changeTarget/diameter", this);
    fChangeTargetDiameterCmd -> SetGuidance("Change the diameter value of the target. "
	                                   "\nDefault value is 7. mm."
					   "\nThe range is between 0 and 15 mm.");
    fChangeTargetDiameterCmd -> SetParameterName("TargetDiameter", true);
    fChangeTargetDiameterCmd -> SetRange("TargetDiameter > 0. && TargetDiameter < 15.");
    fChangeTargetDiameterCmd -> SetDefaultValue(7.*mm);
    fChangeTargetDiameterCmd -> AvailableForStates(G4State_Idle);
    fChangeTargetDiameterCmd -> SetDefaultUnit("mm");
    fChangeTargetDiameterCmd -> SetUnitCandidates("mm");

    // Change Target parameters
    fChangeTargetMaterial = new G4UIdirectory("/changeTarget/designedMaterial/");
    fChangeTargetMaterial -> SetGuidance("Change the Target material choosing isotopes and elements, and their abundance in the target");
    
    //Change target material defining isotopes
    fTargetIsotopeName = new G4UIcmdWithAString("/changeTarget/designedMaterial/isotopeName",this);
    fTargetIsotopeName->SetGuidance("name of the isotope - ex : Ni64");
    fTargetIsotopeName->SetParameterName("IsotopeName",false);
    fTargetIsotopeName->AvailableForStates(G4State_Idle);
    
    fTargetIsotopeZ = new G4UIcmdWithADouble("/changeTarget/designedMaterial/isotopeZ",this);
    fTargetIsotopeZ-> SetGuidance("Z of the isotope");
    fTargetIsotopeZ->SetParameterName("IsotopeZ",false);
    fTargetIsotopeZ->AvailableForStates(G4State_Idle);

    fTargetIsotopeN = new G4UIcmdWithAnInteger("/changeTarget/designedMaterial/isotopeN",this);
    fTargetIsotopeN->SetGuidance("N (number of nucleons) of the isotope");
    fTargetIsotopeN->SetParameterName("IsotopeN",false);
    fTargetIsotopeN->AvailableForStates(G4State_Idle);

    fTargetIsotopeA = new G4UIcmdWithADouble("/changeTarget/designedMaterial/isotopeA",this);
    fTargetIsotopeA->SetGuidance("A of the isotope, in g/cm3");
    fTargetIsotopeA->SetParameterName("IsotopeA",false);
    fTargetIsotopeA->AvailableForStates(G4State_Idle);
    
    //Define elements
    fTargetElementName= new G4UIcmdWithAString("/changeTarget/designedMaterial/ElementName",this);
    fTargetElementName->SetGuidance("Name of the material - ex : PureNi64");
    fTargetElementName->SetParameterName("ElementName",false);
    fTargetElementName->AvailableForStates(G4State_Idle);

    fTargetElementSymbole=new G4UIcmdWithAString("/changeTarget/designedMaterial/ElementSymbole",this);
    fTargetElementSymbole->SetGuidance("Symbole of the element : ex 64Ni");
    fTargetElementSymbole->SetParameterName("ElementSymbole", false);
    fTargetElementSymbole->AvailableForStates(G4State_Idle);

    fTargetElementNComponents = new G4UIcmdWithAnInteger("/changeTarget/designedMaterial/ElementNComponents",this);
    fTargetElementNComponents->SetGuidance("Number of isotopes in the element");
    fTargetElementNComponents->SetParameterName("ElementNComponent", false);
    fTargetElementNComponents->AvailableForStates(G4State_Idle);
    
    fTargetElementAbundance = new G4UIcmdWithADouble("/changeTarget/designedMaterial/IsotopeAbundanceInElement",this);
    fTargetElementAbundance->SetGuidance("Abundance of the isotope in the target");
    fTargetElementAbundance->SetParameterName("IsotopeAbundance",false);
    fTargetElementAbundance->AvailableForStates(G4State_Idle);
    
    //Change material properties
    fChangeTargetMaterialDensityCmd = new G4UIcmdWithADouble("/changeTarget/designedMaterial/materialDensity", this);
    fChangeTargetMaterialDensityCmd -> SetGuidance("Change the density value of the Target Material."
					    "\nDefault value : 8.85 g/cm3.");
    fChangeTargetMaterialDensityCmd -> SetParameterName("TargetMaterialDensity", true);
    fChangeTargetMaterialDensityCmd -> SetDefaultValue(8.85);
    fChangeTargetMaterialDensityCmd -> AvailableForStates(G4State_Idle);
    
    fTargetMaterialNComponents = new G4UIcmdWithAnInteger("/changeTarget/designedMaterial/MaterialNComponents",this);
    fTargetMaterialNComponents->SetGuidance("Number of elements in the target material");
    fTargetMaterialNComponents->SetParameterName("MaterialNComponents",false);
    fTargetMaterialNComponents->AvailableForStates(G4State_PreInit,G4State_Idle);

    fTargetMaterialFractionMass= new G4UIcmdWithADouble("/changeTarget/designedMaterial/MaterialFractionMass",this);
    fTargetMaterialFractionMass->SetGuidance("Fraction mass of the element in the material");
    fTargetMaterialFractionMass->SetParameterName("MaterialFractionMass",false);
    fTargetMaterialFractionMass->AvailableForStates(G4State_Idle);

    fTargetMaterialNaturalElement= new G4UIcmdWithAString("/changeTarget/designedMaterial/naturalElementName",this);
    fTargetMaterialNaturalElement->SetGuidance("Add an element using NIST database");
    fTargetMaterialNaturalElement->SetParameterName("NaturalElement",false);
    fTargetMaterialNaturalElement->AvailableForStates(G4State_Idle);

    fTargetMaterialNaturalElementFractionMass= new G4UIcmdWithADouble("/changeTarget/designedMaterial/naturalElementFractionMass",this);
    fTargetMaterialNaturalElementFractionMass->SetGuidance("Add the fraction mass of the natural element");
    fTargetMaterialNaturalElementFractionMass->SetParameterName("NaturalElementFractionMass",false);
    fTargetMaterialNaturalElementFractionMass->AvailableForStates(G4State_Idle);


    fUpdateMaterial = new G4UIcmdWithoutParameter("/changeTarget/designedMaterial/update",this);
    fUpdateMaterial->SetGuidance("Update the material once its components are defined");
    fUpdateMaterial->AvailableForStates(G4State_Idle);

    //Change material using physics NIST
    fChangeTargetMaterialCmd = new G4UIcmdWithAString("/changeTarget/materialNist", this);
    fChangeTargetMaterialCmd -> SetGuidance("Change the material of your target using the NIST database."
					   "\nTo get the list of the available NIST materials, please select 'TargetMaterial->NistMaterialList'."
					   "\nExample of a NIST material : 'G4_Ni'.");
    fChangeTargetMaterialCmd -> SetParameterName("TargetMaterial",false);
    

    //Change Target thickness
     fChangeTargetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/changeTarget/thickness", this);
     fChangeTargetThicknessCmd -> SetGuidance("Change the thickness value of the Target."
					    "\nDefault value : 0.6 mm.");
    fChangeTargetThicknessCmd -> SetParameterName("TargetThickness", true);
    fChangeTargetThicknessCmd -> SetDefaultValue(0.6*mm);
    fChangeTargetThicknessCmd -> AvailableForStates(G4State_Idle);
    fChangeTargetThicknessCmd -> SetDefaultUnit("mm");
    fChangeTargetThicknessCmd -> SetUnitCandidates("mm");
    

    //////////////////////////
    //Change foil parameters//
    //////////////////////////

    fChangeFoil = new G4UIdirectory("/changeFoil/");
    fChangeFoil -> SetGuidance("Change the Foil thickness");

    // Change Foil Thickness
    fChangeFoilThicknessCmd = new G4UIcmdWithADoubleAndUnit("/changeFoil/thickness", this);
    fChangeFoilThicknessCmd -> SetGuidance("Change the thickness value of the foil "
					  "\nThe default value is 0.32 mm.");
    fChangeFoilThicknessCmd -> SetParameterName("FoilThickness", true);
    fChangeFoilThicknessCmd -> SetDefaultValue(.32*mm);
    fChangeFoilThicknessCmd -> AvailableForStates(G4State_Idle);
    fChangeFoilThicknessCmd -> SetDefaultUnit("mm");
    fChangeFoilThicknessCmd -> SetUnitCandidates("mm");
    
    // Change Target material
    fChangeFoilMaterial = new G4UIdirectory("/changeFoil/designedMaterial/");
    fChangeFoilMaterial -> SetGuidance("Change the Foil material choosing isotopes and elements, and their abundance in the foil");
    
    //Change target material defining isotopes
    fFoilIsotopeName = new G4UIcmdWithAString("/changeFoil/designedMaterial/isotopeName",this);
    fFoilIsotopeName->SetGuidance("name of the isotope - ex : Ni64");
    fFoilIsotopeName->SetParameterName("foilIsotopeName",false);
    fFoilIsotopeName->AvailableForStates(G4State_Idle);
    
    fFoilIsotopeZ = new G4UIcmdWithADouble("/changeFoil/designedMaterial/isotopeZ",this);
    fFoilIsotopeZ-> SetGuidance("Z of the isotope");
    fFoilIsotopeZ->SetParameterName("foilIsotopeZ",false);
    fFoilIsotopeZ->AvailableForStates(G4State_Idle);

    fFoilIsotopeN = new G4UIcmdWithAnInteger("/changeFoil/designedMaterial/isotopeN",this);
    fFoilIsotopeN->SetGuidance("N (number of nucleons) of the isotope");
    fFoilIsotopeN->SetParameterName("foilIsotopeN",false);
    fFoilIsotopeN->AvailableForStates(G4State_Idle);

    fFoilIsotopeA = new G4UIcmdWithADouble("/changeFoil/designedMaterial/isotopeA",this);
    fFoilIsotopeA->SetGuidance("A of the isotope, in g/cm3");
    fFoilIsotopeA->SetParameterName("foilIsotopeA",false);
    fFoilIsotopeA->AvailableForStates(G4State_Idle);
    
    //Define elements

    fFoilElementName= new G4UIcmdWithAString("/changeFoil/designedMaterial/ElementName",this);
    fFoilElementName->SetGuidance("Name of the material - ex : PureNi64");
    fFoilElementName->SetParameterName("foilElementName",false);
    fFoilElementName->AvailableForStates(G4State_Idle);

    fFoilElementSymbole=new G4UIcmdWithAString("/changeFoil/designedMaterial/ElementSymbole",this);
    fFoilElementSymbole->SetGuidance("Symbole of the element : ex 64Ni");
    fFoilElementSymbole->SetParameterName("foilElementSymbole", false);
    fFoilElementSymbole->AvailableForStates(G4State_Idle);

    fFoilElementNComponents = new G4UIcmdWithAnInteger("/changeFoil/designedMaterial/ElementNComponents",this);
    fFoilElementNComponents->SetGuidance("Number of isotopes in the element");
    fFoilElementNComponents->SetParameterName("foilElementNComponent", false);
    fFoilElementNComponents->AvailableForStates(G4State_Idle);
    
    fFoilElementAbundance = new G4UIcmdWithADouble("/changeFoil/designedMaterial/IsotopeAbundanceInElement",this);
    fFoilElementAbundance->SetGuidance("Abundance of the isotope in the foil");
    fFoilElementAbundance->SetParameterName("foilIsotopeAbundance",false);
    fFoilElementAbundance->AvailableForStates(G4State_Idle);
    
    //Change material properties
    fChangeFoilMaterialDensityCmd = new G4UIcmdWithADouble("/changeFoil/designedMaterial/materialDensity", this);
    fChangeFoilMaterialDensityCmd -> SetGuidance("Change the density value of the Target Material");
    fChangeFoilMaterialDensityCmd -> SetParameterName("FoilMaterialDensity", true);
    fChangeFoilMaterialDensityCmd -> AvailableForStates(G4State_Idle);
    
    fFoilMaterialNComponents = new G4UIcmdWithAnInteger("/changeFoil/designedMaterial/MaterialNComponents",this);
    fFoilMaterialNComponents->SetGuidance("Number of elements in the target material");
    fFoilMaterialNComponents->SetParameterName("foilMaterialNComponents",false);
    fFoilMaterialNComponents->AvailableForStates(G4State_Idle);

    fFoilMaterialFractionMass= new G4UIcmdWithADouble("/changeFoil/designedMaterial/MaterialFractionMass",this);
    fFoilMaterialFractionMass->SetGuidance("Fraction mass of the element in the material");
    fFoilMaterialFractionMass->SetParameterName("foilMaterialFractionMass",false);
    fFoilMaterialFractionMass->AvailableForStates(G4State_Idle);

    fFoilMaterialNaturalElement= new G4UIcmdWithAString("/changeFoil/designedMaterial/naturalElementName",this);
    fFoilMaterialNaturalElement->SetGuidance("Add an element using NIST database");
    fFoilMaterialNaturalElement->SetParameterName("foilNaturalElement",false);
    fFoilMaterialNaturalElement->AvailableForStates(G4State_Idle);

    fFoilMaterialNaturalElementFractionMass= new G4UIcmdWithADouble("/changeFoil/designedMaterial/naturalElementFractionMass",this);
    fFoilMaterialNaturalElementFractionMass->SetGuidance("Add the fraction mass of the natural element");
    fFoilMaterialNaturalElementFractionMass->SetParameterName("foilNaturalElementFractionMass",false);
    fFoilMaterialNaturalElementFractionMass->AvailableForStates(G4State_Idle);


    fUpdateFoilMaterial = new G4UIcmdWithoutParameter("/changeFoil/designedMaterial/update",this);
    fUpdateFoilMaterial->SetGuidance("Update the material once its components are defined");
    fUpdateFoilMaterial->AvailableForStates(G4State_Idle);

    //Change foil material using physics NIST

    fChangeFoilMaterialCmd = new G4UIcmdWithAString("/changeFoil/materialNist", this);
    fChangeFoilMaterialCmd -> SetGuidance("Change the material of your foil using the NIST database."
					 "\nTo get the list of the available NIST materials, please select 'TargetMaterial->NistMaterialList'."
					 "\nExample of a NIST material : 'G4_Ni'.");
    fChangeFoilMaterialCmd -> SetParameterName("FoilMaterial",false);
    
   }

STCyclotronDetectorMessenger::~STCyclotronDetectorMessenger()
{
    delete fChangeTarget; 
    delete fChangeTargetDiameterCmd; 
    delete fChangeTargetMaterial;
    delete fTargetIsotopeName;
    delete fTargetIsotopeZ ;
    delete fTargetIsotopeN;
    delete fTargetIsotopeA;
    delete fTargetElementName;
    delete fTargetElementSymbole;
    delete fTargetElementNComponents;
    delete fTargetElementAbundance ;
    delete fChangeTargetMaterialDensityCmd ;
    delete fTargetMaterialNComponents;
    delete fTargetMaterialFractionMass;
    delete fTargetMaterialNaturalElement;
    delete fTargetMaterialNaturalElementFractionMass;
    delete fUpdateMaterial;
    delete fChangeTargetMaterialCmd;
    delete fChangeFoilMaterial;
    delete fFoilIsotopeName;
    delete fFoilIsotopeZ ;
    delete fFoilIsotopeN;
    delete fFoilIsotopeA;
    delete fFoilElementName;
    delete fFoilElementSymbole;
    delete fFoilElementNComponents;
    delete fFoilElementAbundance ;
    delete fChangeFoilMaterialDensityCmd ;
    delete fFoilMaterialNComponents;
    delete fFoilMaterialFractionMass;
    delete fFoilMaterialNaturalElement;
    delete fFoilMaterialNaturalElementFractionMass;
    delete fUpdateFoilMaterial;
    delete fChangeFoilMaterialCmd;
    delete fChangeTargetThicknessCmd;
    delete fChangeFoil; 
    delete fChangeFoilThicknessCmd;
    
}

void STCyclotronDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  //TARGET
  //DIAMETER
  if( command == fChangeTargetDiameterCmd)
    {
      G4double updatedValue  = fChangeTargetDiameterCmd -> GetNewDoubleValue(newValue);
      fDet -> SetTargetDiameter(updatedValue);
    }
 
  //MATERIAL
  else if(command == fTargetIsotopeName)
    {
      fDet -> SetTargetIsotopeName(newValue);
    }

  else if(command == fTargetIsotopeZ)
    {
      fDet -> SetTargetIsotopeZ(fTargetIsotopeZ->GetNewDoubleValue(newValue));
    }

  else if(command == fTargetIsotopeN)
    {
      fDet -> SetTargetIsotopeN(fTargetIsotopeN->GetNewIntValue(newValue));
    }

  else if(command == fTargetIsotopeA)
    {
      fDet -> SetTargetIsotopeA(fTargetIsotopeA->GetNewDoubleValue(newValue));
    }

  else if(command == fTargetElementName)
    {
      fDet -> SetTargetElementName(newValue);
    }

  else if(command == fTargetElementSymbole)
    {
      fDet -> SetTargetElementSymbole(newValue);
    }

  else if(command == fTargetElementNComponents)
    {
      fDet -> SetTargetElementNComponents(fTargetElementNComponents->GetNewIntValue(newValue));
    }

  else if(command == fTargetElementAbundance)
    {
      fDet -> SetTargetElementAbundance(fTargetElementAbundance->GetNewDoubleValue(newValue));
    }


  else if (command == fChangeTargetMaterialDensityCmd )
    {
      G4double updatedValue = fChangeTargetMaterialDensityCmd -> GetNewDoubleValue(newValue);
      fDet -> SetTargetMaterialDensity(updatedValue);
    }

  else if(command == fTargetMaterialNComponents)
    {
      fDet -> SetTargetMaterialNComponents(fTargetMaterialNComponents->GetNewIntValue(newValue));
    }

  else if(command == fTargetMaterialFractionMass)
    {
      fDet -> SetTargetMaterialFractionMass(fTargetMaterialFractionMass->GetNewDoubleValue(newValue));
    }

  else if(command == fUpdateMaterial)
    {
      fDet -> UpdateMaterial();
    }

           //NATURAL ELEMENT
  else if(command == fTargetMaterialNaturalElement)
    {
      fDet ->SetTargetNaturalElement(newValue);
    }

  else if(command == fTargetMaterialNaturalElementFractionMass)
    {
      fDet ->SetTargetNaturalMaterialFractionMass(fTargetMaterialNaturalElementFractionMass->GetNewDoubleValue(newValue));
    }


         //NATURAL MATERIAL
  
  else if (command == fChangeTargetMaterialCmd )
    {
      fDet -> SetTargetMaterial(newValue);
    }

      //THICKNESS

  else if (command == fChangeTargetThicknessCmd )
    {
      G4double updatedValue = fChangeTargetThicknessCmd -> GetNewDoubleValue(newValue);
      fDet -> SetTargetThickness(updatedValue);
    }

  //FOIL

  else if (command == fChangeFoilThicknessCmd )
    {
      G4double updatedValue = fChangeFoilThicknessCmd -> GetNewDoubleValue(newValue);
      fDet -> SetFoilThickness(updatedValue);
    }

   //MATERIAL FOIL
  else if(command == fFoilIsotopeName)
    {
      fDet -> SetFoilIsotopeName(newValue);
    }

  else if(command == fFoilIsotopeZ)
    {
      fDet -> SetFoilIsotopeZ(fFoilIsotopeZ->GetNewDoubleValue(newValue));
    }

  else if(command == fFoilIsotopeN)
    {
      fDet -> SetFoilIsotopeN(fFoilIsotopeN->GetNewIntValue(newValue));
    }

  else if(command == fFoilIsotopeA)
    {
      fDet -> SetFoilIsotopeA(fFoilIsotopeA->GetNewDoubleValue(newValue));
    }

  else if(command == fFoilElementName)
    {
      fDet -> SetFoilElementName(newValue);
    }

  else if(command == fFoilElementSymbole)
    {
      fDet -> SetFoilElementSymbole(newValue);
    }

  else if(command == fFoilElementNComponents)
    {
      fDet -> SetFoilElementNComponents(fFoilElementNComponents->GetNewIntValue(newValue));
    }

  else if(command == fFoilElementAbundance)
    {
      fDet -> SetFoilElementAbundance(fFoilElementAbundance->GetNewDoubleValue(newValue));
    }


  else if (command == fChangeFoilMaterialDensityCmd )
    {
      G4double updatedValue = fChangeFoilMaterialDensityCmd -> GetNewDoubleValue(newValue);
      fDet -> SetFoilMaterialDensity(updatedValue);
    }

  else if(command == fFoilMaterialNComponents)
    {
      fDet -> SetFoilMaterialNComponents(fFoilMaterialNComponents->GetNewIntValue(newValue));
    }

  else if(command == fFoilMaterialFractionMass)
    {
      fDet -> SetFoilMaterialFractionMass(fFoilMaterialFractionMass->GetNewDoubleValue(newValue));
    }

  else if(command == fUpdateFoilMaterial)
    {
      fDet -> UpdateFoilMaterial();
    }

  //NATURAL ELEMENT
  else if(command == fFoilMaterialNaturalElement)
    {
      fDet ->SetFoilNaturalElement(newValue);
    }

  else if(command == fFoilMaterialNaturalElementFractionMass)
    {
      fDet ->SetFoilNaturalMaterialFractionMass(fFoilMaterialNaturalElementFractionMass->GetNewDoubleValue(newValue));
    }

 //NATURAL MATERIAL 
  else if (command == fChangeFoilMaterialCmd )
    {
      fDet -> SetFoilMaterial(newValue);
    }
}
