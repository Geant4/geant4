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
#ifndef MATERIALS_HH
#define MATERIALS_HH

#include "globals.hh"

class G4Material;


class Materials {

 public:
   virtual ~Materials();
   static Materials* Instance();
   virtual void Destroy();
   
   G4Material* GetMaterial(G4String matName);
   
 protected:
   // virtual ~Materials();
   Materials();

   Materials(const Materials& only);
   const Materials& operator=(const Materials& only);

 private:
   static Materials* instance;

   G4Material* hydrogen;
   G4Material* beryllium;
   G4Material* graphite; 
   G4Material* magnesium;
   G4Material* aluminium;
   G4Material* silicon;
   G4Material* liquidArgon;  
   G4Material* titanium;
   G4Material* iron; 
   G4Material* cobalt;
   G4Material* nickel;
   G4Material* indium;
   G4Material* tin; 
   G4Material* copper; 
   G4Material* zinc;  
   G4Material* gallium;
   G4Material* germanium;
   G4Material* zirconium;
   G4Material* molybdenium;
   G4Material* silver;
   G4Material* cadmium;
   G4Material* cesium; 
   G4Material* samarium;
   G4Material* ytterbium; 
   G4Material* tantalum;
   G4Material* tungsten;
   G4Material* gold; 
   G4Material* lead;
   G4Material* uranium;
   G4Material* water; 
   G4Material* quartz; 
   G4Material* ossigeno;
   G4Material* air; 
   G4Material* vacuum;
   G4Material* nytrogen;
};

#endif // MATERIALS_HH
