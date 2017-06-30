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
// $Id: G4FilecoutDestination.hh 103582 2017-04-18 17:24:45Z adotti $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// Implements a cout destination to a file.

//      ---------------- G4FilecoutDestination ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4FILECOUTDESTINATION_HH_
#define G4FILECOUTDESTINATION_HH_

#include <fstream>
#include <memory>

#include "G4coutDestination.hh"

class G4FilecoutDestination : public G4coutDestination
{
  public:

    explicit G4FilecoutDestination(const G4String& fname ,
             std::ios_base::openmode mode = std::ios_base::app )
      : m_name(fname), m_mode(mode), m_output(nullptr) {}
    virtual ~G4FilecoutDestination();

    void SetFileName(const G4String& fname) { m_name = fname; }

    void Open(std::ios_base::openmode mode=std::ios_base::app);
      // By default append to existing file
    void Close();

    virtual G4int ReceiveG4cout(const G4String& msg) override;
    virtual G4int ReceiveG4cerr(const G4String& msg) override;

  private:

    G4String m_name;
    std::ios_base::openmode m_mode;
    std::unique_ptr<std::ofstream> m_output;
};

#endif /* G4FILECOUTDESTINATION_HH_ */
