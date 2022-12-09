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
// G4tgrVolumeAssembly
//
// Class description:
//
// Class for keeping the information of assembly volumes.

// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------
#ifndef G4tgrVolumeAssembly_hh
#define G4tgrVolumeAssembly_hh 1

#include "globals.hh"
#include "G4tgrVolume.hh"
#include "G4ThreeVector.hh"

class G4tgrVolumeAssembly : public G4tgrVolume
{
  public:

    G4tgrVolumeAssembly();
    G4tgrVolumeAssembly(const std::vector<G4String>& wl);
    ~G4tgrVolumeAssembly();

    virtual G4tgrPlace* AddPlace(const std::vector<G4String>& wl);
      // Add a position with the data read from a ':place_assembly' tag

    const G4String& GetComponentName(G4int ii) const
    {
      return theComponentNames[ii];
    }
    const G4String& GetComponentRM(G4int ii) const {return theComponentRMs[ii];}
    G4ThreeVector GetComponentPos(G4int ii) const {return theComponentPos[ii];}
    G4int GetNoComponents() const { return (G4int)theComponentNames.size(); }

    friend std::ostream& operator<<(std::ostream& os,
                                    const G4tgrVolumeAssembly& obj);

  protected:

    std::vector<G4String> theComponentNames;
    std::vector<G4String> theComponentRMs;
    std::vector<G4ThreeVector> theComponentPos;
};

#endif
