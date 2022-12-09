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

// Base class for common implementation for histograms managers.
//
// Author: Ivana Hrivnacova, 10/08/2022  (ivana@ipno.in2p3.fr)

#ifndef G4VTBaseHnManager_h
#define G4VTBaseHnManager_h 1

#include "G4HnInformation.hh"
#include "globals.hh"

#include <vector>
#include <memory>
#include <array>

class G4HnManager;

template <unsigned int DIM>
class G4VTBaseHnManager
{
  // Disable using the object managers outside
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;

  public:
    G4VTBaseHnManager() = default;
    virtual ~G4VTBaseHnManager() = default;

    // deleted copy constructor & assignment operator
    G4VTBaseHnManager(const G4VTBaseHnManager& rhs) = delete;
    G4VTBaseHnManager& operator=(const G4VTBaseHnManager& rhs) = delete;

    // Methods for handling histograms
    virtual G4int Create(const G4String& name, const G4String& title,
               const std::array<G4HnDimension, DIM>& bins,
               const std::array<G4HnDimensionInformation, DIM>& hnInfo) = 0;

    // virtual G4int Create(const G4String& name, const G4String& title,
    //            const std::array<std::vector<G4double>, DIM> edges,
    //            const std::array<G4HnDimensionInformation, DIM>& hnInfo) = 0;

    virtual G4bool Set(G4int id,
               const std::array<G4HnDimension, DIM>& bins,
               const std::array<G4HnDimensionInformation, DIM>& hnInfo) = 0;

    // virtual G4bool Set(G4int id,
    //            const std::array<std::vector<G4double>, DIM>& edges,
    //            const std::array<G4HnDimensionInformation, DIM>& hnInfo) = 0;

    virtual G4bool Scale(G4int id, G4double factor) = 0;

    // Methods to fill histograms
    virtual G4bool Fill(G4int id, std::array<G4double, DIM> value, G4double weight = 1.0) = 0;

    // Access methods
    virtual G4int  GetId(const G4String& name, G4bool warn = true) const = 0;

    // Access to bins parameters
    virtual G4int    GetNbins(unsigned int idim, G4int id) const = 0;
    virtual G4double GetMinValue(unsigned int idim, G4int id) const = 0;
    virtual G4double GetMaxValue(unsigned int idim, G4int id) const = 0;
    virtual G4double GetWidth(unsigned int idim, G4int id) const = 0;

    // Setters for attributes for plotting
    virtual G4bool SetTitle(G4int id, const G4String& title) = 0;
    virtual G4bool SetAxisTitle(unsigned int idim, G4int id, const G4String& title) = 0;

    // Access attributes for plotting
    virtual G4String GetTitle(G4int id) const = 0;
    virtual G4String GetAxisTitle(unsigned int idim, G4int id) const = 0;

    // Methods to list/print histograms
    virtual G4bool WriteOnAscii(std::ofstream& output) = 0;
    virtual G4bool List(std::ostream& output, G4bool onlyIfActive = true) = 0;

    virtual std::shared_ptr<G4HnManager> GetHnManager() = 0;
    virtual const std::shared_ptr<G4HnManager> GetHnManager() const = 0;
};

#endif
