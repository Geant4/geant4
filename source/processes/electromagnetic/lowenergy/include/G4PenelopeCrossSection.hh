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
// $Id: G4PenelopeCrossSection.hh 74626 2013-10-17 07:00:59Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 18 Mar 2010   L. Pandola   1st implementation. 
// 09 Mar 2012   L. Pandola   Add public method (and machinery) to return 
//                            the absolute and the normalized shell cross 
//                            sections independently.
//
// -------------------------------------------------------------------
//
// Class description:
// This class is a container for cross sections and transport momenta 
// calculated by Penelope models (ionisation, bremsstrahlung). It stores 
// PhysicsTables/PhysicsVectors of
// a) the "hard quantities" (above the threshold), 0-th order (cross section) 
//     1-st order (= stopping XS), 2-nd order (= straggling XS)
// b) the "soft quantities" (below threshold), 0-th order (cross section) 
//     1-st order (= stopping XS), 2-nd order (= straggling XS)
// c) total hard cross sections for individual oscillators
// vs. energy. Two versions are available, one with normalized values 
// (good for sampling) and one with absolute values.
// 
// The interface *always* uses energy and cross sections, while internally 
// log(energy) and log(XS) are used.
//
// One instance per each cut-material couple should be created by the 
// calling class. 
//
// Public method to retrieve hard cross section, soft stopping power, 
// total cross section and hard shell cross sections.
//
// Notice: all quantities stored here are *per molecule*
//
// -------------------------------------------------------------------

#ifndef G4PENELOPECROSSSECTION_HH
#define G4PENELOPECROSSSECTION_HH 1

#include "globals.hh"

class G4PhysicsTable;
class G4DataVector;

class G4PenelopeCrossSection
{

public:
  //constructor: one has to give the number of points in each PhysicsVector 
  //(= dimension of the energy grid) and the number of shells (0 is the 
  //default).
  G4PenelopeCrossSection(size_t nOfEnergyPoints,size_t nOfShells=0);
  //
  ~G4PenelopeCrossSection();

  //! Returns total cross section at the given energy
  G4double GetTotalCrossSection(G4double energy) const;
  //! Returns hard cross section at the given energy
  G4double GetHardCrossSection(G4double energy) const;
  //! Returns the total stopping power due to soft collisions
  G4double GetSoftStoppingPower(G4double energy) const;
  //! Returns the hard cross section for the given shell (per molecule)
  G4double GetShellCrossSection(size_t shellID,G4double energy) const;
  //! Returns the hard cross section for the given shell (normalized to 1)
  G4double GetNormalizedShellCrossSection(size_t shellID,G4double energy) const;

  size_t GetNumberOfShells() const {return numberOfShells;};

  //!
  //! Public interface for the master thread
  //! 
  void AddCrossSectionPoint(size_t binNumber, 
                            G4double energy,G4double XH0, G4double XH1, 
			    G4double XH2,
			    G4double XS0, G4double XS1, G4double XS2);
  void AddShellCrossSectionPoint(size_t binNumber,
			         size_t shellID,G4double energy,G4double xs);
  void NormalizeShellCrossSections();

private:
  G4PenelopeCrossSection & operator=(const G4PenelopeCrossSection &right);
  G4PenelopeCrossSection(const G4PenelopeCrossSection&);
  
  G4bool isNormalized;

  size_t numberOfEnergyPoints;
  size_t numberOfShells;

  //all tables are log. XS vs. log E

  //XS0, XS1, XS2 in Penelope nomenclature
  G4PhysicsTable* softCrossSections;

  //XH0, XH1, XH2 in Penelope nomenclature
  G4PhysicsTable* hardCrossSections;
  
  //XS for individual shells
  G4PhysicsTable* shellCrossSections;
  G4PhysicsTable* shellNormalizedCrossSections;

};

#endif

