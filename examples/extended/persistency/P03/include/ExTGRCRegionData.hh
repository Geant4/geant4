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
/// \file persistency/P03/include/ExTGRCRegionData.hh
/// \brief Definition of the ExTGRCRegionData class
//

#ifndef ExTGRCRegionData_h
#define ExTGRCRegionData_h

#include <vector>
#include "globals.hh"

/// Stores cuts per region data
///
/// Changes:     creation   May 2007
/// \author      P. Arce

class ExTGRCRegionData 
{ 
  public:

    ExTGRCRegionData( const std::vector<G4String>& rd );
   ~ExTGRCRegionData();

    void SetCutsData( const std::vector<G4String>& cu );

    //--- Get methods

    G4String GetRegionName() const     { return fRegionName;  }
    std::vector<G4String> GetLVNames() { return fLVNames;     }
    G4double GetGammaCut() const       { return fGammaCut;    }
    G4double GetElectronCut() const    { return fElectronCut; }
    G4double GetPositronCut() const    { return fPositronCut; }
    G4bool CutsAreSet() const          { return fbCutsSet;       }

  private:

    G4String fRegionName;
    std::vector<G4String> fLVNames;
    G4double fGammaCut, fElectronCut, fPositronCut;
    G4bool fbCutsSet;
};

#endif
