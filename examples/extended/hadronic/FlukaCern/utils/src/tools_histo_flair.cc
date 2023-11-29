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
///  \file tools_histo_flair.cc
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

#include "tools_histo_flair.hh"

#include <fstream>
#include <iomanip>

#include <stdlib.h>
#include <stdarg.h>
#include <cstring>
#include <memory>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace tools {
  namespace histo {
    namespace flair {

      // ***************************************************************************
      // Returns abscissa name.
      // ***************************************************************************
      G4String getAbscissaName(const Abscissa abscissaKind) {
        G4String abscissaName = "";
        if (abscissaKind == Abscissa::KineticEnergy) {
          abscissaName = "kE";
        }
        else if (abscissaKind == Abscissa::Z) {
          abscissaName = "Z";
        }
        else if (abscissaKind == Abscissa::A) {
          abscissaName = "A";
        }
        return abscissaName;
      }


      // ***************************************************************************
      // Returns abscissa unit.
      // ***************************************************************************
      G4String getAbscissaUnit(const Abscissa abscissaKind) {
        G4String abscissaUnit = "";
        if (abscissaKind == Abscissa::KineticEnergy) {
          abscissaUnit = "GeV";
        }
        return abscissaUnit;
      }
      

      // ***************************************************************************
      // Returns abscissa unit value.
      // ***************************************************************************
      G4double getAbscissaUnitValue(const Abscissa abscissaKind) {
        G4double abscissaUnitValue = 1.;
        if (abscissaKind == Abscissa::KineticEnergy) {
          abscissaUnitValue = G4UnitDefinition::GetValueOf(getAbscissaUnit(Abscissa::KineticEnergy));
        }
        return abscissaUnitValue;
      }


      // ***************************************************************************
      // String format like printf() for strings.
      // This sformat function is extracted as-is from geoviewer, 
      // and was written by V. Vlachoudis in 2010.
      // ***************************************************************************
      G4String sformat(const char* fmt, ...) {
        std::unique_ptr<char[]> formatted;
        G4int n = (G4int)strlen(fmt) * 2; // Reserve two times as much as the length of the fmt_str.
        while (true) {
          // Wrap the plain char array into the unique_ptr.
          formatted.reset(new char[n]);
          strcpy(&formatted[0], fmt);
          va_list ap;
          va_start(ap, fmt);
          G4int final_n = vsnprintf(&formatted[0], n, fmt, ap);
          va_end(ap);
          if (final_n < 0 || final_n >= n) {
            n += std::abs(final_n - n + 1);
          }
          else {
            break;
          }
        }
        return G4String(formatted.get());
      }
    

      // ***************************************************************************
      // Returns MC error.
      // ***************************************************************************
      G4double computeError(const G4double mean, 
                          const G4double sumSquares, 
                          const G4int numEvents) {

        const G4double var  = std::max(0., sumSquares / numEvents - mean*mean);
        const G4double err  = ((mean != 0. && numEvents > 1) ? 
                             std::sqrt(var / (numEvents-1)) / mean 
                             : 1.);

        return err * 100.;
      }


      // ***************************************************************************
      // Dump G4H1 (histo mode) to a Flair-compatible format.
      // ***************************************************************************
      void dumpG4H1HistoInFlairFormat(std::ofstream& output,
                                      const G4int indexInOutputFile,
                                      const G4String& histoName,
                                      G4H1* const histo,
                                      const Abscissa abscissaKind,
                                      const G4String& binSchemeName,
                                      const G4int numEvents,
                                      const G4double sumSquaredEventTotals,
                                      const G4double sumSquaredEventInRangeTotals) {

        const G4bool isProfile = false;
        dumpG4H1InFlairFormat(output,
                              indexInOutputFile,
                              histoName,
                              histo,
                              abscissaKind,
                              binSchemeName,
                              numEvents,
                              isProfile,
                              sumSquaredEventTotals,
                              sumSquaredEventInRangeTotals);
      }


      // ***************************************************************************
      // Dump G4H1 (profile mode, ie no stats) to a Flair-compatible format.
      // ***************************************************************************
      void dumpG4H1ProfileInFlairFormat(std::ofstream& output,
                                        const G4int indexInOutputFile,
                                        const G4String& histoName,
                                        G4H1* const histo,
                                        const Abscissa abscissaKind,
                                        const G4String& binSchemeName) {

        const G4bool isProfile = true;
        const G4int numEvents = 1;
        dumpG4H1InFlairFormat(output,
                              indexInOutputFile,
                              histoName,
                              histo,
                              abscissaKind,
                              binSchemeName,
                              numEvents,
                              isProfile);
      }


      // ***************************************************************************
      // Dump G4H1 (profile + histo modes) to a Flair-compatible format.
      // ***************************************************************************
      void dumpG4H1InFlairFormat(std::ofstream& output,
                                 const G4int indexInOutputFile,
                                 const G4String& histoName,
                                 G4H1* const histo,
                                 const Abscissa abscissaKind,
                                 const G4String& binSchemeName,
                                 const G4int numEvents,
                                 const G4bool isProfile,
                                 const G4double sumSquaredEventTotals,
                                 const G4double sumSquaredEventInRangeTotals) {

        const G4int numBins = histo->axis().bins();
        const G4double minAbscissa = histo->axis().lower_edge();
        const G4double maxAbscissa = histo->axis().upper_edge();

        const G4String abscissaName = getAbscissaName(abscissaKind);
        const G4String abscissaUnit = getAbscissaUnit(abscissaKind);
        const G4double abscissaUnitValue = getAbscissaUnitValue(abscissaKind);


        // START HEADER.
        if (indexInOutputFile != 1) { output << G4endl << G4endl; }
        output << std::left << std::setw(13) << "# Detector: " << indexInOutputFile 
               << " " << histoName << G4endl;
        output << std::left << std::setw(13) << "# Dim: " << histo->get_dimension() << G4endl;
        output << std::left << std::setw(13) << "# Entries: " << numEvents << G4endl;

        if (!isProfile) {
          // Total integral
          const G4double total = static_cast<G4double>(histo->sum_all_bin_heights()) / numEvents;   
          const G4double totalError = computeError(total, 
                                                 sumSquaredEventTotals, 
                                                 numEvents);
          output << std::left << std::setw(13) << "# Total: "
                 << sformat("%-15g [1/pr] +/- %6.2f [%] \n",
                                       total,
                                       totalError);

          // In range integral
          const G4double inRange = static_cast<G4double>(histo->sum_bin_heights()) / numEvents;
          const G4double inRangeError = computeError(inRange, 
                                                   sumSquaredEventInRangeTotals, 
                                                   numEvents);
          output << std::left << std::setw(13) << "# InRange: "
                 << sformat("%-15g [1/pr] +/- %6.2f [%] \n",
                                       inRange,
                                       inRangeError);

          // Underflow
          const G4double underflow = static_cast<G4double>(histo->bins_sum_w().at(0)) / numEvents;
          const G4double underflowError = computeError(underflow, 
                                                     histo->bins_sum_w2().at(0), 
                                                     numEvents);
          output << std::left << std::setw(13) << "# Under: "
                 << sformat("%-15g [1/pr] +/- %6.2f [%] \n",
                                       underflow,
                                       underflowError);

          // Overflow
          const G4double overflow = static_cast<G4double>(histo->bins_sum_w().at(numBins+1)) / numEvents;
          const G4double overflowError = computeError(overflow, 
                                                    histo->bins_sum_w2().at(numBins+1), 
                                                    numEvents);
          output << std::left << std::setw(13) << "# Over: "
                 << sformat("%-15g [1/pr] +/- %6.2f [%] \n",
                                       overflow,
                                       overflowError);
        }

        // END HEADER.
        output << std::left << std::setw(13) << "# Log: " << binSchemeName << G4endl;
        output << std::left << std::setw(13) << "# Histogram: "
               << numBins
               << " " << minAbscissa / abscissaUnitValue
               << " " << maxAbscissa / abscissaUnitValue
               << "  [" << abscissaUnit << "]"
               << G4endl;

        output << std::left << std::setw(13) << "# Column: " 
               << "1  " << abscissaName << "_min  " 
               << "[" << abscissaUnit << "]" 
               << G4endl;
        output << std::left << std::setw(13) << "# Column: " 
               << "2  " << abscissaName << "_max  "
               << "[" << abscissaUnit << "]" 
               << G4endl;
        output << std::left << std::setw(13) << "# Column: "
               << "3  dN/d(" << abscissaName << ")  "
               << "[1" << (abscissaUnit.empty() ? "" : ("/" + abscissaUnit)) << "/pr]" 
               << G4endl;
        if (!isProfile) {
          output << std::left << std::setw(13) << "# Column: " 
                 << "4  error  "
                 << "[%]"
                 << G4endl;
        }

    
        // HISTO / PROFILE DATA.
        for (G4int binIndex = 0; binIndex < numBins; ++binIndex) {

          const G4double mean = histo->bin_height(binIndex) / numEvents;
          const G4double err = computeError(mean, 
                                          histo->bins_sum_w2().at(binIndex+1), 
                                          numEvents);

          // histo mode
          if (!isProfile) {
            output << sformat("%15g %15g %15g %6.2f\n",
                                         histo->axis().bin_lower_edge(binIndex) / abscissaUnitValue,
                                         histo->axis().bin_upper_edge(binIndex) / abscissaUnitValue,
                                         mean / histo->axis().bin_width(binIndex) * abscissaUnitValue,
                                         err);
          }
          // profile mode
          else {
            output << sformat("%15g %15g %15g\n",
                                         histo->axis().bin_lower_edge(binIndex) / abscissaUnitValue,
                                         histo->axis().bin_upper_edge(binIndex) / abscissaUnitValue,
                                         mean / histo->axis().bin_width(binIndex) * abscissaUnitValue);
          }
        } // loop on all bins
      } // dumpG4H1InFlairFormat


    } // namespace flair
  } // namespace histo
} // namespace tools

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
