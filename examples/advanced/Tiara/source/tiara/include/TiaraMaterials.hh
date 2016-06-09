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
// $Id: TiaraMaterials.hh,v 1.4 2003/06/25 09:12:46 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
//
// Class TiaraMaterials
//

#ifndef TiaraMaterials_hh
#define TiaraMaterials_hh TiaraMaterials_hh

#include "globals.hh"
#include <map>

class G4Material;
class G4Element;

typedef std::map< G4String , G4Element* > TiaraMapSymbolElement;
typedef std::map< G4String, G4Material* > TiaraMapNameMaterial;

class TiaraMaterials{
public:
  TiaraMaterials();
  ~TiaraMaterials();

  
  G4Material *GetMaterial(const G4String &matName) const;
  
  G4Material *CreateAir();
  G4Material *CreateConcrete();
  G4Material *CreateMCNPConcrete();
  G4Material *CreateIron();
  G4Material *CreateVakuum();

private:
  TiaraMaterials(const TiaraMaterials &);

  void FillElementMap(const G4String &name, 
		      const G4String &symbol,
		      G4int Z,
		      G4double A);

  TiaraMaterials &operator=(const TiaraMaterials &);

  TiaraMapSymbolElement fMapSymbolElement;
  TiaraMapNameMaterial fMapNameMaterial;
};

#endif
