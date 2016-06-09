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
// $Id: TiaraVComponent.hh,v 1.4 2003/06/25 09:12:53 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraVComponent
//

#ifndef TiaraVComponent_hh
#define TiaraVComponent_hh TiaraVComponent_hh

#include <vector>
#include "G4ThreeVector.hh"
#include "TiaraPart.hh"

class G4LogicalVolume;

typedef std::vector<TiaraPart> TiaraParts;

class TiaraVComponent {
public:
  TiaraVComponent();
  virtual ~TiaraVComponent();

  virtual TiaraParts GetParts() = 0;
};

#endif
