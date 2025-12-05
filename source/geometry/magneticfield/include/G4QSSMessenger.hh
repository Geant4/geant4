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
// G4QSSMessenger
//
// Messenger for QSS Integrator driver

// Author: Leandro Gomez Vidal (Univ. Buenos Aires), October 2021
// --------------------------------------------------------------------
#ifndef G4QSSMessenger_HH
#define G4QSSMessenger_HH

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

#include "G4QSSParameters.hh"

class G4QSSMessenger : public G4UImessenger
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4QSSMessenger();
    ~G4QSSMessenger() override;

    /**
     * Applies command to the associated object.
     */
    void SetNewValue(G4UIcommand* command, G4String newValues) override;

    /* Hacky - much easier to access G4QSSMessenger from G4QSSDriver than the other way around
     * Multithreading seems to cause weird stuff with Driver/Stepper instances, maybe making
     * some thread-local copies or something, outside of construction? */

    static G4QSSMessenger* instance();

    /**
     * Accessors.
     */
    inline G4int GetQssOrder() { return G4QSSParameters::Instance()->GetQssOrder(); }
    inline G4double Get_dQRel() { return G4QSSParameters::Instance()->Get_dQRel(); }
    inline G4double Get_dQMin() { return G4QSSParameters::Instance()->Get_dQMin(); }
    inline G4int GetMaxSubsteps() { return G4QSSParameters::Instance()->GetMaxSubsteps(); }   

    enum StepperSelection
    {
      None = 1,
      G4QSS2 = 2,
      G4QSS3 = 3,
      NumMethods,
    };

    /**
     * Stepper selection, G4QSS2 or G4QSS3.
     */
    void selectStepper(const std::string&);
    StepperSelection selectedStepper();

    /**
     * Sets QSS order. To be suppressed in favour of the method
     * it calls in G4QSSParameters.
     */
    G4bool SetQssOrder(G4int order);
   
  private:

    /**
     * Internal methods -- could be suppressed in future.
     */
    G4bool Set_dQMin( G4double dvalue );
    G4bool Set_dQRel( G4double value );
    G4bool SetMaxSubsteps( G4int number );
   
  private:

    StepperSelection _selectedStepper= StepperSelection::None;

    G4UIdirectory*             qssCmdDir;
    G4UIcmdWithADoubleAndUnit* dQMinCmd;
    G4UIcmdWithADouble*        dQRelCmd;
    G4UIcmdWithAString*        stepperSelectorCmd;
    G4UIcmdWithAnInteger*      maxSubstepsCmd;
};

#endif
