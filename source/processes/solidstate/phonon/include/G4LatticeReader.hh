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
/// \file processes/phonon/include/G4LatticeReader.hh
/// \brief Definition of the G4LatticeReader class
//
// NOTE:  This reader class for logical lattices should be moved to
//	  materials/ after the 10.0 release (and this comment removed).
// $Id: G4LatticeReader.hh 76799 2013-11-15 20:30:53Z mkelsey $
//
// 20131115  Move ctor, dtor implementations to .cc file.

#ifndef G4LatticeReader_h
#define G4LatticeReader_h 1

#include "globals.hh"
#include <iosfwd>

class G4LatticeLogical;

class G4LatticeReader {
public:
  G4LatticeReader(G4int vb=0);
  ~G4LatticeReader();

  void SetVerboseLevel(G4int vb) { verboseLevel = vb; }

  G4LatticeLogical* MakeLattice(const G4String& filepath);

protected:
  G4bool OpenFile(const G4String& filepath);
  G4bool ProcessToken();
  G4bool ProcessValue(const G4String& name);	// Numerical parameters
  G4bool ProcessConstants();			// Four dynamical constants
  G4bool ProcessMap();				// Velocity magnitudes file
  G4bool ProcessNMap();				// Direction vectors file
  G4bool ReadMapInfo();				// Get map file parameters
  G4bool SkipComments();				// Everything after '#'
  void CloseFile();

private:
  G4int verboseLevel;		// For reporting progress, also use G4VERBOSE

  std::ifstream* psLatfile;	// Configuration file being read
  G4LatticeLogical* pLattice;	// Lattice under construction (not owned)

  G4String fMapPath;		// Path to config file to find velocity maps
  G4String fToken;		// Reusable buffers for reading file
  G4double fValue;		// ... floating point data values
  G4String fMap, fsPol;		// ... map filename and polarization code
  G4int    fPol, fNX, fNY;	// ... map binning in each direction

  static const G4String fDataDir;	// Directory path ($G4LATTICEDATA)
};

#endif	/* G4LatticeReader_h */
