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
// Author: M.A. Cortes-Giraldo, Universidad de Sevilla
//
// History changelog prior creation of this example:
// - 17/10/2009: inheritance removed (not needed)
// - 17/10/2009: version 1.0
// - 07/10/2024: version 2.0 for Geant4 example (MT compliant)
// - 18/10/2025: version 3.0
//   - Creation of IAEASourceIdRegistry for thread-safe source_id assignation
//

#ifndef G4IAEAphspWriter_hh
#define G4IAEAphspWriter_hh 1

#include "G4ThreeVector.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4Run;
class G4Step;

class G4IAEAphspWriterStack;

//------------------------------------------------------------------------------

class G4IAEAphspWriter
{

public:

  G4IAEAphspWriter(const G4String filename);
  ~G4IAEAphspWriter();

  void OpenIAEAphspOutFiles(const G4Run*);
  void WriteIAEAParticle(const size_t idx, const G4int nStat, const G4int pdg,
			 const G4double kinE, const G4double wt,
			 const G4ThreeVector pos, const G4ThreeVector momDir);
  void CloseIAEAphspOutFiles();

  void AddZphsp(const G4double zphsp);
  void SetDataFromIAEAStack(const G4IAEAphspWriterStack* );
  // void UpdateHeaders();

  void SetFileName(const G4String name)     { fFileName = name; }
  void SetConstVariable(G4int idx, G4double value);
  void SumOrigHistories(size_t idx, G4int value)
  { fOrigHistories->at(idx) += value; }

  const G4String GetFileName() const                 { return fFileName; }
  const std::vector<G4double>* GetZphspVec() const   { return fZphspVec; }
  const std::vector<G4int>* GetOrigHistoriesVec() const
  { return fOrigHistories; }


private:

  G4IAEAphspWriter() = default;
  void WriteIAEAParticle(const G4Step* aStep, const G4int zStopIdx);

  // ------------
  // DATA MEMBERS
  // ------------

  // FILE PROPERTIES

  G4String fFileName;
  // Must include the path but not any of the IAEA extensions.

  std::vector<G4double>* fZphspVec = nullptr;
  // Vector storing the z-value of the phsp planes.

  std::map<G4int, G4double>* fConstVariables = nullptr;
  // Map to store the value of the variables set as constant for the phsp's
  // to save disk space.

  // COUNTERS AND FLAGS

  std::vector<G4int>* fOrigHistories = nullptr;
  // Vector bookkeeping the number of original histories recorded for each phsp.
  // For regular simulations, all the elements should have the same value.

  G4int fIAEASourcesOpen;
  // This is a counter to consider if there is another IAEA file already open,
  // e.g. a source IAEAphsp file to generate particles.

};

#endif
