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
/// \file radiobiology/include/Manager.hh
/// \brief Definition of the RadioBio::Manager class

#ifndef RadiobiologyManager_h
#define RadiobiologyManager_h 1

#include "VRadiobiologicalQuantity.hh"

#include <map>
namespace RadioBio
{

class ManagerMessenger;

class Manager
{
  private:
    Manager() = default;

  public:
    // Singleton structure
    ~Manager();
    static Manager* CreateInstance();
    static Manager* GetInstance();

    // Initialization
    void InitializeAll();

    // Calculation
    void ComputeAll();

    // Clear accumulated data
    void ResetAll();

    // Output to text files (called at the end of run)
    void StoreAll();

    // Print manager parameters
    void PrintParameters();

    // Add from the proper accumulable
    void DigestAccumulables();

    // Register a Quantity
    bool Register(VRadiobiologicalQuantity*, G4String);

    // Hook for a quantity given the name
    VRadiobiologicalQuantity* GetQuantity(G4String);

  private:
    // Map with all the radiobiological quantities.
    // The string is the key to access the proper class
    std::map<G4String, VRadiobiologicalQuantity*> fQuantities = {};

    inline static Manager* fInstance = nullptr;

  
};

}  // namespace RadioBio

#endif  // RadioBioManager_h
