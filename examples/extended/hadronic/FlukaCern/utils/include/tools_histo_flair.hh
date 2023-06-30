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
///  \file tools_histo_flair.hh
///  \brief Tools to dump any G4H1 into Flair-compatible format.
//
//  Author: G.Hugo, 08 December 2022
//
// ***************************************************************************
//
//      tools_histo_flair
//
///  Tools to dump any G4H1 into Flair-compatible format.
///
///  These tools are fully application-agnostic, 
///  and could also be placed outside of the G4 examples, 
///  as a core G4 Analysis Manager extension.
//
// ***************************************************************************

#ifndef TOOLS_HISTO_FLAIR_HH
#define TOOLS_HISTO_FLAIR_HH

#include "globals.hh"
#include "g4hntools_defs.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace tools {
  namespace histo {
    namespace flair {

      // Supported abscissa types (can be extended)
      enum Abscissa {
        KineticEnergy,
        Z,
        A
      };
      
      G4String getAbscissaName(const Abscissa abscissaKind);
      G4String getAbscissaUnit(const Abscissa abscissaKind);
      G4double getAbscissaUnitValue(const Abscissa abscissaKind);

      G4String sformat(const char* fmt, ...);

      G4double computeError(const G4double mean, 
                          const G4double sumSquares, 
                          const G4int numEvents);

      // Dump G4H1 (histo mode) to a Flair-compatible format.
      void dumpG4H1HistoInFlairFormat(std::ofstream& output,
                                      const G4int indexInOutputFile,
                                      const G4String& histoName,
                                      G4H1* const histo,
                                      const Abscissa abscissaKind,
                                      const G4String& binSchemeName,
                                      const G4int numEvents,
                                      const G4double sumSquaredEventTotals = 0.,
                                      const G4double sumSquaredEventInRangeTotals = 0.);

      // Dump G4H1 (profile mode, ie no stats) to a Flair-compatible format.
      void dumpG4H1ProfileInFlairFormat(std::ofstream& output,
                                        const G4int indexInOutputFile,
                                        const G4String& histoName,
                                        G4H1* const histo,
                                        const Abscissa abscissaKind,
                                        const G4String& binSchemeName);

      void dumpG4H1InFlairFormat(std::ofstream& output,
                                 const G4int indexInOutputFile,
                                 const G4String& histoName,
                                 G4H1* const histo,
                                 const Abscissa abscissaKind,
                                 const G4String& binSchemeName,
                                 const G4int numEvents,
                                 const G4bool isProfile = false,
                                 const G4double sumSquaredEventTotals = 0.,
                                 const G4double sumSquaredEventInRangeTotals = 0.);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
