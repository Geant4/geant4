//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VVolumeMaterialScanner.hh,v 1.2 2005/12/13 08:43:56 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// class G4VVolumeMaterialScanner
//
// Class description:
//
// Interface for repeated volumes or parameterisations that are able to 
//  tabulate their materials.

// History:
//  18.11.05  J.Apostolakis  First version

#ifndef G4VVOLUMEMATERIALSCANNER_HH
#define G4VVOLUMEMATERIALSCANNER_HH
class G4Material;
// #include "G4Material.hh"

class G4VVolumeMaterialScanner
{
  public:
    G4VVolumeMaterialScanner() {}
    virtual ~G4VVolumeMaterialScanner() {}

    virtual G4int       GetNumberOfMaterials() const =0;
    virtual G4Material* GetMaterial(G4int idx) const =0;
};
   
#endif
 
