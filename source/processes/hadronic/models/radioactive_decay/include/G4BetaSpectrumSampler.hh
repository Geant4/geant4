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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4BetaSpectrumSampler.hh                                          //
//  Author: D.H. Wright                                                       //
//  Date:   21 November 2022                                                  //
//  Description: samples a spectrum which is a piece-wise linear function of  //
//               energy. The CDF is calculated by trapezoidal integration     //
//               and within bins the line y = mx + b is sampled.              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4BetaSpectrumSampler_h
#define G4BetaSpectrumSampler_h 1

#include "globals.hh"
#include "Randomize.hh"
#include <vector>

class G4BetaSpectrumSampler
{
  public:
    G4BetaSpectrumSampler(const G4double* aPDF, G4int pdfSize, G4double e);

    ~G4BetaSpectrumSampler() = default;

    G4double shoot();

  private:
    G4double sampleSlopedLine();  // Method to sample under sloped line
  
    std::vector<double> pdf;      // Spectrum shape
    std::vector<double> cdf;      // Cumulative distribution function

    G4double eEnd;                // Endpoint energy of spectrum
    G4int nBins;                  // Number of lower bin edges in spectrum 
    G4int lowerBinEdge;           // Index of lower bin edge 
    G4int upperBinEdge;           // Index of upper bin edge
    G4double ylower{0.0};         // Value at lower bin edge
    G4double yupper{0.0};         // Value at upper bin edge
};
#endif

