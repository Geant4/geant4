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

#ifndef G4IRTUtils_hh
#define G4IRTUtils_hh

#include "globals.hh"
#include <memory>

class G4ErrorFunction;
class G4IRTUtils {

public:
	G4IRTUtils() = default;
   ~G4IRTUtils() = default;
    static G4double EffectiveDistance(const G4double& rc,
                                      const G4double& r0);

    static G4double GetKact(const G4double& obs, 
                            const G4double& dif)
    {
        return (obs == 0 || dif == 0) ? 0 : dif * obs/(dif - obs);
    }
    
    static G4double GetRCutOff();
    static G4double GetRCutOff(G4double tCutOff);
    static G4double GetDNADistanceCutOff();
};



#endif

