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
// $Id: G4ToolsAnalysisReader.hh 91012 2015-06-15 10:38:07Z ihrivnac $

// The common implementation of analysis reader classes.

// Author: Ivana Hrivnacova, 20/07/2015  (ivana@ipno.in2p3.fr)

#ifndef G4ToolsAnalysisReader_h
#define G4ToolsAnalysisReader_h 1

#include "G4VAnalysisReader.hh"
#include "globals.hh"

#include "tools/histo/h1d" 
#include "tools/histo/h2d" 
#include "tools/histo/h3d" 
#include "tools/histo/p1d" 
#include "tools/histo/p2d" 
#include "tools/wroot/ntuple"

class G4H1ToolsManager;
class G4H2ToolsManager;
class G4H3ToolsManager;
class G4P1ToolsManager;
class G4P2ToolsManager;

namespace tools {
namespace histo {
class hmpi;  
}  
}
 
class G4ToolsAnalysisReader : public  G4VAnalysisReader
{
  public:
    explicit G4ToolsAnalysisReader(const G4String& type, G4bool isMaster = true);
    virtual ~G4ToolsAnalysisReader();
    
    // Access methods
    tools::histo::h1d*  GetH1(G4int id, G4bool warn = true) const;
    tools::histo::h2d*  GetH2(G4int id, G4bool warn = true) const;
    tools::histo::h3d*  GetH3(G4int id, G4bool warn = true) const;
    tools::histo::p1d*  GetP1(G4int id, G4bool warn = true) const;
    tools::histo::p2d*  GetP2(G4int id, G4bool warn = true) const;
    
  protected:
    // methods
    G4bool Reset();

    // data members
    G4H1ToolsManager*  fH1Manager;
    G4H2ToolsManager*  fH2Manager;
    G4H3ToolsManager*  fH3Manager;
    G4P1ToolsManager*  fP1Manager;
    G4P2ToolsManager*  fP2Manager;
 };

#include "G4ToolsAnalysisReader.icc"

#endif

