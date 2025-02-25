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
// G4MSSteppingAction implementation
//
// Author: M.Asai, 5 May 2006
// --------------------------------------------------------------------

#include "G4MSSteppingAction.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4Region.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"

#include <map>

// --------------------------------------------------------------------
void G4MSSteppingAction::Initialize(G4bool rSens, G4Region* reg)
{
  regionSensitive = rSens;
  theRegion = reg;
  length = 0.;
  x0 = 0.;
  lambda = 0.;
  shape_mat_info_v.clear();
}

// --------------------------------------------------------------------
void G4MSSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Region* region = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();

  if (regionSensitive && (region != theRegion)) return;

  G4double stlen = aStep->GetStepLength();
  const G4Material* material = preStepPoint->GetMaterial();
  length += stlen;
  x0 += stlen / (material->GetRadlen());
  lambda += stlen / (material->GetNuclearInterLength());

  // store information per step (1 geantino step = 1 solid)
  {
    shape_mat_info_v.push_back({});

    shape_mat_info_t & thisMaterialInfo = shape_mat_info_v.back();

    // calculate average mass number and atomic number
    {
      const std::vector<const G4Element*> * ElementVector_ptr = material->GetElementVector();
      for( auto & element : *ElementVector_ptr)
      {
        thisMaterialInfo.aveA += element->GetA();
        thisMaterialInfo.aveZ += element->GetZ();
      }
      thisMaterialInfo.aveA /= ElementVector_ptr->size();
      thisMaterialInfo.aveZ /= ElementVector_ptr->size();
    }
    thisMaterialInfo.density            = material->GetDensity();
    thisMaterialInfo.radiation_length   = material->GetRadlen();
    thisMaterialInfo.interaction_length = material->GetNuclearInterLength();
    thisMaterialInfo.thickness            = aStep->GetStepLength();
    thisMaterialInfo.integrated_thickness = length;
    thisMaterialInfo.lambda            = stlen / (material->GetNuclearInterLength());
    thisMaterialInfo.x0                = stlen / (material->GetRadlen());
    thisMaterialInfo.integrated_lambda = lambda;
    thisMaterialInfo.integrated_x0     = x0;
    thisMaterialInfo.entry_point = preStepPoint->GetPosition();
    thisMaterialInfo.exit_point  = aStep->GetPostStepPoint()->GetPosition();
    thisMaterialInfo.material_name = material->GetName();
  }

}

void G4MSSteppingAction::PrintEachMaterialVerbose(std::ostream & oss)
{

  G4int colwidth = 11;
  G4int matname_colwidth = 15;

  oss << G4endl<< G4endl;
  oss << "   Material         Atomic properties (averaged)   Radiation Interaction             Integr.    Lambda      X0                          Entry point                              Exit point\n     name            Mass        Z       density    length     length    Thickness  Thickness                                               (cm)                                     (cm)\n                  (g/mole)               (g/cm3)     (cm)       (cm)       (cm)        (cm)       \n";
  oss << G4endl;


  std::ios::fmtflags os_flags (oss.flags());
  for( auto & matInfo : shape_mat_info_v)
  {
     oss << std::setw(matname_colwidth) << std::left << matInfo.GetName(matname_colwidth) << "   ";

     oss << std::setw(colwidth) << std::fixed << std::setprecision(2) << matInfo.aveA / (CLHEP::g/CLHEP::mole);
     oss << std::setw(colwidth) << std::fixed << std::setprecision(2) << matInfo.aveZ;
     oss << std::setw(colwidth) << std::scientific << std::setprecision(2) << matInfo.density / (CLHEP::g/CLHEP::cm3);
     oss << std::setw(colwidth) << std::scientific << matInfo.radiation_length / CLHEP::cm;
     oss << std::setw(colwidth) << std::scientific << matInfo.interaction_length / CLHEP::cm;
     oss << std::setw(colwidth) << std::scientific << matInfo.thickness/CLHEP::cm;
     oss << std::setw(colwidth) << std::scientific << matInfo.integrated_thickness/CLHEP::cm;
     oss << std::setw(colwidth) << std::scientific << matInfo.lambda;
     oss << std::setw(colwidth) << std::scientific << matInfo.x0;
     oss << std::setw(colwidth) << " ";
     oss << std::scientific <<  std::right << "(";
     oss << std::scientific <<  std::left  << matInfo.entry_point.x()/CLHEP::cm;
     oss <<  ", "  << matInfo.entry_point.y()/CLHEP::cm;
     oss <<  ", "  << matInfo.entry_point.z()/CLHEP::cm << ")";
     oss << std::setw(colwidth) << " ";
     oss << std::scientific <<  std::right << "(";
     oss << std::scientific <<  std::left  << matInfo.exit_point.x()/CLHEP::cm;
     oss << ", "  << matInfo.exit_point.y()/CLHEP::cm;
     oss << ", "  << matInfo.exit_point.z()/CLHEP::cm << ")";
     oss << G4endl;
     oss << G4endl;
  }
  oss.flags(os_flags);  // Restore original stream format
}

void G4MSSteppingAction::PrintIntegratedMaterialVerbose(std::ostream& oss)
{
  // create database (db) of material name (key) and information
  std::map<G4String, shape_mat_info_t> mat_db;
  // accumulate information for each material name into mat_db
  for(auto & matInfo : shape_mat_info_v)
  {
     G4String key = matInfo.material_name;
     if( 0 == mat_db.count( key ) )
     {
       mat_db[key] = shape_mat_info_t{};
     }

    mat_db[key].x0 += matInfo.x0;
    mat_db[key].thickness += matInfo.thickness;
    mat_db[key].lambda += matInfo.lambda;
  }

  std::ios::fmtflags os_flags (oss.flags());
  oss << std::scientific << std::setprecision(2) << '\t';
  for(auto & [key,mat] : mat_db)
    oss << '\t' << key
        << '\t'<< mat.thickness/CLHEP::mm
        << '\t'<< mat.x0
        << '\t'<< mat.lambda;
  oss.flags(os_flags);  // Restore original stream format
  return;
}
