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
// $Id: G4tgrFileReader.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrFileReader
//
// Class description:
//
// This service provides access to detector description data.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrFileReader_H
#define G4tgrFileReader_H 1

#include "globals.hh"

#include <vector>

class G4tgrVolume;
class G4tgrVolumeMgr;
class G4tgrLineProcessor;

class G4tgrFileReader
{
  public:  // with description

    static G4tgrFileReader* GetInstance();  
      // Get the only instance 

    virtual ~G4tgrFileReader();

    G4bool ReadFiles();

    void AddTextFile(const G4String& fname) { theTextFiles.push_back( fname ); }
    void SetLineProcessor( G4tgrLineProcessor* lp ) { theLineProcessor = lp; }
    G4tgrLineProcessor* GetLineProcessor() const { return theLineProcessor; }

  protected:

    G4tgrFileReader();

  private:

    static G4ThreadLocal G4tgrFileReader* theInstance;  

    std::vector<G4String> theTextFiles;
    G4tgrLineProcessor* theLineProcessor;

};

# endif 
