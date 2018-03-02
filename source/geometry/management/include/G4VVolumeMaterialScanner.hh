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
// $Id: G4VVolumeMaterialScanner.hh 108480 2018-02-15 14:28:27Z gcosmo $
//
// class G4VVolumeMaterialScanner
//
// Class description:
//
// Interface for repeated volumes or parameterisations that are able to 
//  tabulate their materials.

// History:
//  18.11.05 J.Apostolakis Initial version
// --------------------------------------------------------------------
#ifndef G4VVOLUMEMATERIALSCANNER_HH
#define G4VVOLUMEMATERIALSCANNER_HH

#include "G4Types.hh"

class G4Material;

class G4VVolumeMaterialScanner
{
  public:
    G4VVolumeMaterialScanner() {;}
    virtual ~G4VVolumeMaterialScanner() {;}

    virtual G4int       GetNumberOfMaterials() const = 0;
    virtual G4Material* GetMaterial(G4int idx) const = 0;
};
   
#endif
