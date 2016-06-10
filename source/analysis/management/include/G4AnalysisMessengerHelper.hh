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
// $Id: G4AnalysisMessengerHelper.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The helper class for histogram and profiles messengers.
// It implements reusable commands in /analysis/[hn|pn] directory.
//
// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisMessengerHelper_h
#define G4AnalysisMessengerHelper_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class G4VAnalysisManager;
class G4UIdirectory;
class G4UIcommand;

class G4AnalysisMessengerHelper
{
  public:
    // types
    struct BinData {
      BinData() :     
        fNbins(0),
        fVmin(0.),
        fVmax(0.),
        fSunit(""),
        fSfcn(""),
        fSbinScheme("") {}
      G4int    fNbins;
      G4double fVmin;
      G4double fVmax;
      G4String fSunit;
      G4String fSfcn;
      G4String fSbinScheme;
    };
    // types
    struct ValueData {
      ValueData() :     
        fVmin(0.),
        fVmax(0.),
        fSunit(""),
        fSfcn("") {}
      G4double fVmin;
      G4double fVmax;
      G4String fSunit;
      G4String fSfcn;
    };

  public:
    // Make available utility Update method 
    friend class G4HnMessenger;

 public:
    explicit G4AnalysisMessengerHelper(const G4String& hnType);
    ~G4AnalysisMessengerHelper();
   
    // methods to create commands
    std::unique_ptr<G4UIdirectory>  CreateHnDirectory() const; 

    std::unique_ptr<G4UIcommand>  CreateSetTitleCommand(
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetBinsCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetValuesCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;
    std::unique_ptr<G4UIcommand>  CreateSetAxisCommand(const G4String& axis,
                                        G4UImessenger* messenger) const;

    // methods to read command paremeters
    void GetBinData(BinData& data, std::vector<G4String>& parameters, 
                    G4int& counter) const;
    void GetValueData(ValueData& data, std::vector<G4String>& parameters, 
                    G4int& counter) const;

    // warnings
    void WarnAboutParameters(G4UIcommand* command, G4int nofParameters) const;
    void WarnAboutSetCommands() const;

  private:
    // methods
    G4String Update(const G4String& str, const G4String& axis = "") const;

    // data members
    G4String  fHnType;
};
  
#endif

