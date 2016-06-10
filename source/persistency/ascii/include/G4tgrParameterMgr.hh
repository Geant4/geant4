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
// $Id: G4tgrParameterMgr.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrParameterMgr
//
// Class description:
//
// Class to manage parameters. It is a singleton.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrParameterMgr_h
#define G4tgrParameterMgr_h

#include "globals.hh"

#include <map>
#include <vector>

typedef std::map< G4String, G4String > G4mapss;

class G4tgrParameterMgr 
{ 
  public:   // with description

    static G4tgrParameterMgr* GetInstance();  
      // Get the only instance 

    void AddParameterNumber( const std::vector<G4String>& wl,
                                   G4bool mustBeNew = 0 );
    void AddParameterString( const std::vector<G4String>& wl,
                                   G4bool mustBeNew = 0 );
      // Add to theParameterList

    void CheckIfNewParameter( const std::vector<G4String>& wl,
                                    G4bool mustBeNew  );
      // Check if it is new and that there are 3 words

    G4String FindParameter( const G4String& name, G4bool exists = true );
     // Find a Parameter with name 'name'. 

    void DumpList();
      // Dump list of parameters

  private:

    G4tgrParameterMgr();
   ~G4tgrParameterMgr();

  private:

    G4mapss theParameterList;
      // Map of Parameter's: G4String is the Parameter name,
      // double is its value

    static G4ThreadLocal G4tgrParameterMgr* theInstance;
};

#endif
