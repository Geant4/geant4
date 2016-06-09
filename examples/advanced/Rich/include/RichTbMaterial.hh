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
// Rich advanced example for Geant4
// RichTbMaterial.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////

#ifndef RichTbMaterial_h
#define RichTbMaterial_h 1
#include "globals.hh"
#include <vector>
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4OpticalSurface.hh"
#include "RichTbRunConfig.hh"
#include "AerogelTypeSpec.hh"

class RichTbMaterial {
  
public:  
  RichTbMaterial();

  RichTbMaterial(RichTbRunConfig*);
  virtual ~RichTbMaterial();
  G4double ConvertAgelRIndex(G4double,G4int);
  G4Material* getAir() {return RichTbAmbientAir;}
  G4Material* getTAir() {return RichTbTubeAir;}
  G4Material* getNitrogenGas() {return RichTbNitrogenGas;}
  G4Material* getH2O() {return RichTbH2O;}
  G4Material* getMirrorQuartz() {return RichTbMirrorQuartz;}
  G4Material* getQuartzWMaterial() {return  RichTbQuartzWindowMaterial;}
  G4Material* getHpdQuartzWMaterial() {return  HpdQuartzWindowMaterial;}
  G4Material* getPadHpdQuartzWMaterial() {return  PadHpdQuartzWindowMaterial;}
  G4Material* getRichTbFilterMaterial(G4int FilterTNum )
  {return RichTbFilterMaterial[FilterTNum] ; }
  G4Material* getGlassD263FilterMaterial() {return  GlassD263FilterMaterial;}
  G4Material* getAerogelMaterial(G4int AgelNum ) 
  {return RichTbAerogelMaterial[AgelNum]; }
  G4Material* getAerogelTypeA() {return RichTbAerogelTypeA; }
  G4Material* getAerogelTypeB() {return RichTbAerogelTypeB; }
  G4Material* getAerogelTypeC() {return RichTbAerogelTypeC; }
  G4Material* getAerogelTypeD() {return RichTbAerogelTypeD; }
  G4Material* getAerogelTypeE() {return RichTbAerogelTypeE; }
  G4Material* getC4F10() {return RichTbC4F10; }
  G4Material* getCF4() {return RichTbCF4; }
  G4Material* getVacuum() {return RichTbVacuum;}
  G4Material* getAluminium() {return RichTbAluminium; } 
  G4Material* getHpdTubeMaterial() {return  HpdTubeMaterial; }
  G4Material* getHpdSiDetMaterial() {return HpdSiDetMaterial;}
  G4Material* getHpdSiCoatingMaterial() {return HpdSiCoatingMaterial;}
  G4Material* getPadHpdPhCathodeMaterial() {return PadHpdPhCathodeMaterial;}

  G4Material* getPlasticAg() {return RichTbPlasticAg;}
  G4OpticalSurface* getOpticalMirrorSurface() 
  {return RichTbOpticalMirrorSurface; } 
  G4OpticalSurface* getOpticalEnclosureSurface() 
  {return RichTbOpticalEnclosureSurface; } 
  G4OpticalSurface* getOpticalGasQuartzWSurface() 
  {return RichTbOpticalGasQuartzWSurface; } 
  G4OpticalSurface* getOpticalHpdMetalSurface()
  {return  RichTbOpticalHpdMetalSurface;}
  G4OpticalSurface* getOpticalHpdTQuartzWSurface()
  {return  HpdTQuartzWSurface; }
  G4OpticalSurface* getOpticalHpdQuartzWPhCathodeSurface()
  {return  HpdQuartzWPhCathodeSurface; }
  G4OpticalSurface* getOpticalPhCathodeSkinSurface()
  {return  PhCathodeSkinSurface; }
  G4OpticalSurface* getOpticalPhCathodeBorderSurface()
  {return  PhCathodeBorderSurface; }
  G4OpticalSurface* getOpticalFilterSurface()
  {return  RichTbOpticalFilterSurface; }
  G4OpticalSurface* getHpdSiCoatSurface()
  {return  HpdSiCoatSurface; }

private:
  G4Material* RichTbAmbientAir;
  G4Material* RichTbTubeAir;
  G4Material* RichTbNitrogenGas;
  G4Material* RichTbH2O;
  G4Material* RichTbMirrorQuartz;
  G4Material* RichTbQuartzWindowMaterial;
  std::vector<G4Material*> RichTbAerogelMaterial;

  G4Material* RichTbAerogelTypeA;
  G4Material* RichTbAerogelTypeB;
  G4Material* RichTbAerogelTypeC;
  G4Material* RichTbAerogelTypeD;
  G4Material* RichTbAerogelTypeE;
  G4Material* RichTbC4F10;
  G4Material* RichTbCF4;
  G4Material* RichTbVacuum;
  G4Material* RichTbAluminium;
  G4Material* HpdTubeMaterial;
  G4Material* HpdSiDetMaterial;
  G4Material* HpdSiCoatingMaterial;
  G4Material* HpdQuartzWindowMaterial;
  G4Material* PadHpdQuartzWindowMaterial;

  std::vector<G4Material*> RichTbFilterMaterial;
  G4Material* GlassD263FilterMaterial;
  G4Material* PadHpdPhCathodeMaterial;
  G4Material* RichTbPlasticAg;
  G4OpticalSurface* RichTbOpticalMirrorSurface;
  G4OpticalSurface* RichTbOpticalEnclosureSurface;
  G4OpticalSurface*  RichTbOpticalGasQuartzWSurface;
  G4OpticalSurface* HpdTQuartzWSurface;
  G4OpticalSurface* HpdQuartzWPhCathodeSurface;
  G4OpticalSurface* PhCathodeSkinSurface;
  G4OpticalSurface* PhCathodeBorderSurface;
  G4OpticalSurface* HpdSiCoatSurface;
  //The following metal surface used although it would be redundant.
  // The metal volume has no ref. index.
  G4OpticalSurface* RichTbOpticalHpdMetalSurface;

  G4OpticalSurface* RichTbOpticalFilterSurface;


  RichTbRunConfig* rConfig;
};



#endif 

