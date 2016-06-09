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
// $Id: ExN07GapParameterisation.hh,v 1.1 2003/03/10 01:43:35 asaim Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal size.
//

#ifndef ExN07GapParameterisation_H
#define ExN07GapParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Material;

class ExN07GapParameterisation : public G4VPVParameterisation
{ 
  public:
  
    ExN07GapParameterisation();
    virtual ~ExN07GapParameterisation();
   
    virtual void ComputeTransformation(const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;
    virtual G4Material* ComputeMaterial(const G4int copyNo,
                                G4VPhysicalVolume* physVol);

  private:
    G4int       numberOfLayers;
    G4Material* absMaterial;
    G4Material* gapMaterial;

  public:
    inline void SetNumberOfLayers(G4int nl)
    { numberOfLayers = nl; }
    inline void SetAbsorberMaterial(G4Material* mat)
    { absMaterial = mat; }
    inline void SetGapMaterial(G4Material* mat)
    { gapMaterial = mat; }
};


#endif


