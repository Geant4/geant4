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

// Common implementation for histograms managers.
//
// Author: Ivana Hrivnacova, 10/08/2022  (ivana@ipno.in2p3.fr)

#ifndef G4THnToolsManager_h
#define G4THnToolsManager_h 1

#include "G4VTBaseHnManager.hh"
#include "G4THnManager.hh"
#include "globals.hh"

#include <memory>
#include <array>
#include <string_view>

class G4AnalysisManagerState;
class G4HnManager;
class G4UImessenger;

template <unsigned int DIM, typename HT>
class G4THnToolsManager : public G4VTBaseHnManager<DIM>,
                          public G4THnManager<HT>

{
  // Disable using the object managers outside
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;

  public:
    G4THnToolsManager(const G4AnalysisManagerState& state);
    virtual ~G4THnToolsManager() = default;

    // deleted copy constructor & assignment operator
    G4THnToolsManager(const G4THnToolsManager& rhs) = delete;
    G4THnToolsManager& operator=(const G4THnToolsManager& rhs) = delete;

    // Methods for handling histograms
    G4int Create(const G4String& name, const G4String& title,
               const std::array<G4HnDimension, DIM>& bins,
               const std::array<G4HnDimensionInformation, DIM>& hnInfo) override;

    G4bool Set(G4int id,
               const std::array<G4HnDimension, DIM>& bins,
               const std::array<G4HnDimensionInformation, DIM>& hnInfo) override;

    virtual G4bool Scale(G4int id, G4double factor) override;

    // Methods to fill histograms
    virtual G4bool Fill(G4int id, std::array<G4double, DIM> value, G4double weight = 1.0) override;

    // Access methods
    virtual G4int  GetId(const G4String& name, G4bool warn = true) const  override;

    // Access to bins parameters
    virtual G4int    GetNbins(unsigned int idim, G4int id) const  override;
    virtual G4double GetMinValue(unsigned int idim, G4int id) const override;
    virtual G4double GetMaxValue(unsigned int idim, G4int id) const override;
    virtual G4double GetWidth(unsigned int idim, G4int id) const override;

    // Setters for attributes for plotting
    virtual G4bool SetTitle(G4int id, const G4String& title) override;
    virtual G4bool SetAxisTitle(unsigned int idim, G4int id, const G4String& title)  override;

    // Access attributes for plotting
    virtual G4String GetTitle(G4int id) const  override;
    virtual G4String GetAxisTitle(unsigned int idim, G4int id) const override;

    // Methods to list/print histograms
    virtual G4bool WriteOnAscii(std::ofstream& output) override;
            // Function with specialization per histo type
    virtual G4bool List(std::ostream& output, G4bool onlyIfActive = true) override;

    // // Access to Hn manager
    virtual std::shared_ptr<G4HnManager> GetHnManager() override;
    virtual const std::shared_ptr<G4HnManager> GetHnManager() const override;

  protected:
    // Functions from base class
    using G4THnManager<HT>::RegisterT;
    using G4THnManager<HT>::GetTHnInFunction;
    using G4THnManager<HT>::GetTInFunction;
    using G4THnManager<HT>::GetTId;
    using G4THnManager<HT>::IsVerbose;
    using G4THnManager<HT>::Message;

    static constexpr std::string_view fkClass { "G4THnToolsManager" };

  private:
    void UpdateInformation(G4HnInformation* hnInformation,
                           const std::array<G4HnDimensionInformation, DIM>& hnInfo);

    G4HnInformation* AddInformation(const G4String& name,
                       const std::array<G4HnDimensionInformation, DIM>& hnInfo);

    void AddAnnotation(HT* ht,
                       const std::array<G4HnDimensionInformation, DIM>& hnInfo);

    G4bool CheckName(const G4String& name) const;

    // Functions with specialization per histo type

    HT* CreateToolsHT(const G4String& title,
                      const std::array<G4HnDimension, DIM>& bins,
                      const std::array<G4HnDimensionInformation, DIM>& hnInfo);

    void ConfigureToolsHT(HT* ht,
                      const std::array<G4HnDimension, DIM>& bins,
                      const std::array<G4HnDimensionInformation, DIM>& hnInfo);

    G4bool FillHT(HT* ht, const G4HnInformation& hnInformation,
                     std::array<G4double, DIM>& value, 
                     G4double weight = 1.0);

   // redefine tools::histo::key_axis_x_title() - not available via upper types
   static const std::array<std::string, 3> fkKeyAxisTitle;

   std::unique_ptr<G4UImessenger> fMessenger;
};

// inline functions

#include "G4THnMessenger.hh"
#include "G4THnToolsManager.icc"
#include "G4THnMessenger.icc"
    // include messenger and its implementation here to avoid include recursion

#endif
