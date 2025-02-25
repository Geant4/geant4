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
// G4MSSteppingAction
//
// Class description:
//
// Stepping action for material scanner.

// Author: M.Asai, 5 May 2006
// --------------------------------------------------------------------
#ifndef G4MSSteppingAction_hh
#define G4MSSteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <vector>

class G4Region;


class G4MSSteppingAction : public G4UserSteppingAction
{
  public:
    G4MSSteppingAction() = default;
    ~G4MSSteppingAction() override = default;

    void Initialize(G4bool rSens, G4Region* reg);
    void UserSteppingAction(const G4Step*) override;

    inline G4double GetTotalStepLength() const { return length; }
    inline G4double GetX0() const { return x0; }
    inline G4double GetLambda0() const { return lambda; }

    /// Print material properties verbosely for each step of geantino
    /// This function is useful for single shot scans.
    void PrintEachMaterialVerbose(std::ostream & oss);

    /// Print list of {material name, thickness, x0, lambda}, integrated by material name
    /// This function is useful for global scans.
    void PrintIntegratedMaterialVerbose(std::ostream & oss);

  private:
    G4bool regionSensitive = false;
    G4Region* theRegion = nullptr;
    G4double length = 0.0;
    G4double x0 = 0.0;
    G4double lambda = 0.0;
    struct shape_mat_info_t
    {
      /// Calculated average atomic number
      G4double aveZ               = 0.0;
      /// Calculated average mass number
      G4double aveA               = 0.0;
      /// Density of material given by user
      G4double density            = 0.0;
      /// Material radiation length
      G4double radiation_length   = 0.0;
      /// Material interaction length
      G4double interaction_length = 0.0;
      /// Step of the geantino
      G4double thickness          = 0.0;
      /// Integrated path of the geantino
      G4double integrated_thickness = 0.0;
      /// Calculated x0 = thickness/radiation_length
      G4double x0                 = 0.0;
      /// Calculated lambda = thickness/interaction_length
      G4double lambda             = 0.0;
      /// Integrated x0
      G4double integrated_x0       = 0.0;
      /// Integrated lambda
      G4double integrated_lambda   = 0.0;
      /// Entry point of the geantino into the solid
      G4ThreeVector entry_point   = { };
      /// Exit point of the geantino out of the solid
      G4ThreeVector exit_point    = { };
      /// Material name. Composition is not checked.
      /// That is, if there are two identical materials
      /// except for the name, they are treated as different
      G4String material_name      = { };
      /// Getter that returns the full material name
      const G4String& GetName() { return material_name; }
      /// Getter that returns the material name, splitted in blocks of length 'column_width'
      G4String GetName(G4int column_width)
      {
        auto input_name_length = (G4int)material_name.length();
        if( input_name_length < column_width)  { return material_name; }

        G4String formatted_name;
        for (std::size_t i = 0; i < material_name.length(); i += column_width)
        {
          // for each block of characters of length 'column_width', append '\n'
          formatted_name += material_name.substr(i, column_width);
          if (i + column_width < material_name.length())
          {
            formatted_name += '\n';
          }
          // append spaces for last block of characters so its length corresponds to column_width
          else
          {
            formatted_name+=G4String( column_width-(input_name_length%column_width),' ');
          }
        }
        return formatted_name;
      }
    };
    std::vector<shape_mat_info_t> shape_mat_info_v;

};

#endif
