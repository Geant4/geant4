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
// $Id: ExTGRCRegionData.hh,v 1.3 2010-11-05 08:52:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:      P. Arce
// Changes:     creation   May 2007
// ---------------------------------------------------------------------------

#ifndef ExTGRCRegionData_h
#define ExTGRCRegionData_h

#include <vector>
#include "globals.hh"

class ExTGRCRegionData 
{ 
  public:

    ExTGRCRegionData( const std::vector<G4String>& rd );
   ~ExTGRCRegionData();

    void SetCutsData( const std::vector<G4String>& cu );

    //--- Get methods

    G4String GetRegionName() const     { return theRegionName;  }
    std::vector<G4String> GetLVNames() { return theLVNames;     }
    G4double GetGammaCut() const       { return theGammaCut;    }
    G4double GetElectronCut() const    { return theElectronCut; }
    G4double GetPositronCut() const    { return thePositronCut; }
    G4bool CutsAreSet() const          { return bCutsSet;       }

  private:

    G4String theRegionName;
    std::vector<G4String> theLVNames;
    G4double theGammaCut, theElectronCut, thePositronCut;
    G4bool bCutsSet;
};

#endif
