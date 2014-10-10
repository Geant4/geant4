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
// $Id: G4MTcoutDestination.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4MTcoutDestination.hh
//
// ---------------------------------------------------------------
#ifndef G4MTcoutDestination_H
#define G4MTcoutDestination_H

#include "globals.hh"
#include "G4coutDestination.hh"
#include <iostream>
#include <sstream>
#include <fstream>

class G4MTcoutDestination : public G4coutDestination
{
  public:

    G4MTcoutDestination(const G4int& threadId,
       std::ostream& co=std::cout, std::ostream&  ce=std::cerr);
    virtual ~G4MTcoutDestination();

    virtual G4int ReceiveG4cout(const G4String&);
    virtual G4int ReceiveG4cerr(const G4String&);

    void SetCoutFileName(const G4String& fileN = "G4cout.txt", G4bool ifAppend = true);
    void SetCerrFileName(const G4String& fileN = "G4cerr.txt", G4bool ifAppend = true);
    void EnableBuffering(G4bool flag=true);
    void SetPrefixString(const G4String& wd = "G4WT");
    void SetIgnoreCout(G4int tid = 0);
    void SetIgnoreInit(G4bool val=true) { ignoreInit = val; }
    G4String GetPrefixString() const { return prefix; }
    G4String GetFullPrefixString() const {
        std::stringstream os;
        os<<prefix<<id;
        return os.str();
    }

  private:

    void CloseCoutFile();
    void CloseCerrFile();
    void DumpBuffer();
  
  private:

    std::ostream& finalcout;
    std::ostream& finalcerr;
    const G4int id;
    G4bool useBuffer;
    G4bool threadCoutToFile;
    G4bool threadCerrToFile;
    G4bool ignoreCout;
    G4bool ignoreInit;

    std::ostringstream cout_buffer;
    std::ostringstream cerr_buffer;
    std::ofstream coutFile;
    std::ofstream cerrFile;
    G4String prefix;
};

#endif
