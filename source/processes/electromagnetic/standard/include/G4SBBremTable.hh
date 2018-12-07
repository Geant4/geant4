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
//
// ----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4SBBremTable
//
// Author:        Mihaly Novak
//
// Creation date: 15.07.2018
//
// Modifications:
//
// Class description:
//
// Utility class to handle sampling tables for the Seltzer-Berger scalled brems-
// strahlung differential cross sections. It makes possible fast (significantly
// faster than the rejection) sampling of the emitted photon energy in case of
// interactions. An object from this class is supposed to be a member of the
// Seltzer-Berger model for e-/e+ bremsstrahlung photon emission model. Note,
// that one object from this class can handle both e- and e+ cases (containes
// e+ correction in the SampleEnergyTransfer method only).
//
// ----------------------------------------------------------------------------

#ifndef G4SBBremTable_h
#define G4SBBremTable_h 1

#include "globals.hh"
#include "G4String.hh"

#include <vector>

// forward declar
class G4MaterialCutsCouple;

class G4SBBremTable {

public:
   // CTR/DTR
   G4SBBremTable();

  ~G4SBBremTable();

   // loads and init sampling tables: lowe/highe are the low/high energy usage
   // limits of the corresponding Seltzerberger-model.
   void Initialize(const G4double lowe, const G4double highe);

   // clean away all sampling tables and makes ready for re-initialisation
   void ClearSamplingTables();

   // run-time method to sample energy transferred to the emitted photon
   double SampleEnergyTransfer(const G4double eekin, const G4double leekin,
                               const G4double gcut , const G4double dielSupConst,
                               const G4int    izet , const G4int matCutIndx,
                               const bool     iselectron);

   // used only for development: print out table related information
   // void Dump();

private:

  void  BuildSamplingTables();

  void  InitSamplingTables();

  void  LoadSTGrid();

  void  LoadSamplingTables(G4int iz);

  void  ReadCompressedFile(const G4String &fname, std::istringstream &iss);

private:

  // Sampling-Table point: describes one [E_i],[kappa_j] point
  struct STPoint {
    G4double fCum;    // value of the cumulative function
    G4double fParA;   // rational function approximation based interp. parameter
    G4double fParB;   // rational function approximation based interp. parameter
  };

  // Sampling-Table: describes one [E_j] e- energy point i.e. one Table
  struct STable {
    // cumulative values for the kappa-cuts: kappa_cut_i=E_gamma_cut_i/E_el_j
    std::vector<G4double> fCumCutValues;
    // as many STPoint-s as kappa values
    std::vector<STPoint>  fSTable;
  };

  // Sampling-Tables for a given Z:
  // describes all tables (i.e. for all e- energies) for a given element (Z)
  struct SamplingTablePerZ {
    SamplingTablePerZ() : fNumGammaCuts(0), fMinElEnergyIndx(-1), fMaxElEnergyIndx(-1) {}
    size_t                fNumGammaCuts;     // number of gamma-cut for this
    G4int                 fMinElEnergyIndx;  // max(i) such E_i <= E for all E
    G4int                 fMaxElEnergyIndx;  // min(i) such E_i >= E for all E
    std::vector<STable*>  fTablesPerEnergy;  // as many table as e-ekin grid point
    //the different gamma-cut values that are defined for this element(Z) and ln
    std::vector<G4double> fGammaECuts;
    std::vector<G4double> fLogGammaECuts;
    // the couple index element stores the corresponding (sorted) gamma-cut index
    std::vector<size_t>   fMatCutIndxToGamCutIndx;
    // temporary vector to store some indecis during initialisation
    std::vector< std::vector<size_t> >   fGamCutIndxToMatCutIndx;
  };

  // simple linear search: most of the time faster than anything in our case
  G4int LinSearch(const std::vector<STPoint>& vect,
                  const G4int size,
                  const G4double val);

private:

  // pre-prepared sampling tables are available:
  G4int                           fMaxZet;      // max Z number
  G4int                           fNumElEnergy; // # e- kine (E_k) per Z
  G4int                           fNumKappa;    // # red. photon eners per E_k

  // min/max electron kinetic energy usage limits
  G4double                        fUsedLowEenergy;
  G4double                        fUsedHighEenergy;
  G4double                        fLogMinElEnergy;
  G4double                        fILDeltaElEnergy;

  // e- kinetic energy and reduced photon energy grids and tehir logarithms
  std::vector<G4double>           fElEnergyVect;
  std::vector<G4double>           fLElEnergyVect;
  std::vector<G4double>           fKappaVect;
  std::vector<G4double>           fLKappaVect;

  // container to store samplingtables per Z (size is fMaxZet+1)
  std::vector<SamplingTablePerZ*> fSBSamplingTables;

};

#endif
