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
//

#ifndef G4VCHEMISTRYWORLD_HH
#define G4VCHEMISTRYWORLD_HH

#include <memory>
#include <map>
class G4DNABoundingBox;
class G4Material;
class G4MolecularConfiguration;
class G4VChemistryWorld
{
 public:
    using MolType = const G4MolecularConfiguration*;
    G4VChemistryWorld() = default;
    virtual ~G4VChemistryWorld() = default;

    virtual void ConstructChemistryBoundary() = 0;
    virtual void ConstructChemistryComponents() = 0;

    std::map<MolType,double>::iterator begin()
    {
        return fpChemicalComponent.begin();
    }

    std::map<MolType,double>::iterator end()
    {
      return fpChemicalComponent.end();
    }

    size_t size()
    {
        return fpChemicalComponent.size();
    }

    std::map<MolType,double>::const_iterator begin_const()
    {
        return fpChemicalComponent.begin();
    }

    std::map<MolType,double>::const_iterator end_const()
    {
        return fpChemicalComponent.end();
    }

    G4DNABoundingBox* GetChemistryBoundary() const
    {
        return fpChemistryBoundary.get();
    }
   protected:
    std::unique_ptr<G4DNABoundingBox> fpChemistryBoundary;
    std::map<MolType,double> fpChemicalComponent;
};

#endif
