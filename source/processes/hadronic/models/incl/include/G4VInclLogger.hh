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
// $Id: G4VInclLogger.hh,v 1.2 2010-10-26 02:47:59 kaitanie Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#ifndef G4VInclLogger_hh
#define G4VInclLogger_hh 1

#include "globals.hh"
#include "G4String.hh"

/**
 * INCL logger interface
 *
 * Interface that allows us to generate graphs and histograms of the
 * INCL internal variables.
 */
class G4VInclLogger {

public:
  G4VInclLogger() {};
  ~G4VInclLogger() {};

  /**
   * Book 1D histogram.
   */
  virtual void bookHistogram1D(G4String name, G4int bins, G4double xmin, G4double xmax) = 0;

  /**
   * Book 2D histogram.
   */
  virtual void bookHistogram2D(G4String name, G4int binsx, G4double xmin, G4double xmax,
			       G4int binsy, G4double ymin, G4double ymax) = 0;

  /**
   * Fill 1D histogram.
   */
  virtual void fillHistogram1D(G4String name, G4double value) = 0;

  /**
   * Fill 2D histogram.
   */
  virtual void fillHistogram2D(G4String name, G4double xvalue, G4double yvalue) = 0;

  /**
   * Save the histograms.
   */
  virtual void saveHistograms() = 0;
};

#endif
