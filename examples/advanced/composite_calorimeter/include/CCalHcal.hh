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
///////////////////////////////////////////////////////////////////////////////
// File: CCalHcal.hh
// Description: Equipped to construct the geometry of the hadron calorimeter
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalHcal_h
#define CCalHcal_h 1

#include "CCalDetector.hh"

class CCalHcal: public CCalDetector
{
public:
  //Constructor and Destructor
  CCalHcal(const G4String &name);
  virtual ~CCalHcal();

  //Get Methods
  G4String getGenMat()                      const {return genMaterial;}
  G4double   getDy_2Cal()                   const {return dy_2Cal;}
  G4double   getDx_2Cal()                   const {return dx_2Cal;}
  G4double   getXposCal()                   const {return xposCal;}
  G4String getBoxMat()                      const {return boxMaterial;}
  G4int      getNBox()                      const {return nBox;}
  G4double   getDy_2Box()                   const {return dy_2Box;}
  G4double   getDx_2Box()                   const {return dx_2Box;}
  G4double   getWallThickBox()              const {return wallThickBox;}
  G4double   getXposBox(unsigned int i)     const {return xposBox[i];}
  G4int      getNLayerScnt()                const {return nLayerScnt;}
  G4int      getTypeScnt(unsigned int i)    const {return typeLayerScnt[i];}
  G4int      getMotherScnt(unsigned int i)  const {return mothLayerScnt[i];}
  G4double   getXposScnt(unsigned int i)    const {return xposLayerScnt[i];}
  G4int      getNLayerAbs()                 const {return nLayerAbs;}
  G4int      getTypeAbs(unsigned int i)     const {return typeLayerAbs[i];}
  G4int      getMotherAbs(unsigned int i)   const {return mothLayerAbs[i];}
  G4double   getXposAbs(unsigned int i)     const {return xposLayerAbs[i];}
  G4String getAbsMat()                      const {return absMaterial;}
  G4int      getNAbsorber()                 const {return nAbsorber;}
  G4double   getDy_2Abs(     )              const {return dy_2Absorber;}
  G4double   getDx_2Abs(unsigned int i)     const {return dx_2Absorber[i];}
  G4String getScntMat()                     const {return scntMaterial;}
  G4String getWrapMat()                     const {return wrapMaterial;}
  G4String getPlasMat()                     const {return plasMaterial;}
  G4int      getNScintillator()             const {return nScintillator;}
  G4double   getDy_2ScntLay(unsigned int i) const {return dy_2ScntLayer[i];}
  G4double   getDx_2ScntLay(unsigned int i) const {return dx_2ScntLayer[i];}
  G4double   getDx_2Wrap(unsigned int i)    const {return dx_2Wrapper[i];}
  G4double   getDx_2FrontP(unsigned int i)  const {return dx_2FrontPlastic[i];}
  G4double   getDx_2BackP(unsigned int i)   const {return dx_2BackPlastic[i];}
  G4double   getDx_2Scnt(unsigned int i)    const {return dx_2Scintillator[i];}

protected:
  virtual G4int readFile();
  virtual void constructDaughters();

private:
  G4String genMaterial;              //General material
  G4double   dy_2Cal;                //Half width     of the Hcal 
  G4double   dx_2Cal;                //Half thickness of the Hcal
  G4double   xposCal;                //Position in mother

  G4String boxMaterial;              //Material of boxes
  G4int      nBox;                   //Number of boxes
  G4double   dy_2Box;                //Half width     of the Boxes
  G4double   dx_2Box;                //Half thickness of the Boxes
  G4double   wallThickBox;           //Wall thickness of the boxes
  G4double*  xposBox;                //Position in mother

  G4int      nLayerScnt;             //Number of scintillator layers
  G4int*     typeLayerScnt;          //Layer type
  G4int*     mothLayerScnt;          //Mother type
  G4double*  xposLayerScnt;          //Position in mother

  G4int      nLayerAbs;              //Number of absorber     layers
  G4int*     typeLayerAbs;           //Layer type
  G4int*     mothLayerAbs;           //Mother type
  G4double*  xposLayerAbs;           //Position in mother

  G4String absMaterial;              //Material of absorbers
  G4int      nAbsorber;              //Number of absorber types
  G4double   dy_2Absorber;           //Half width     of the absorbers
  G4double*  dx_2Absorber;           //Half thickness of the absorbers

  G4String scntMaterial;             //Material of Scintillator
  G4String wrapMaterial;             //Material of Wrapper
  G4String plasMaterial;             //Material of plastic cover
  G4int      nScintillator;          //Number of scintillator types
  G4double*  dy_2ScntLayer;          //Half width     of scintillator layers
  G4double*  dx_2ScntLayer;          //Half thickness of scintillator layers
  G4double*  dx_2Wrapper;            //Half thickness of wrappers
  G4double*  dx_2FrontPlastic;       //Half thickness of front plastic
  G4double*  dx_2BackPlastic;        //Half thickness of back  plastic
  G4double*  dx_2Scintillator;       //Half thickness of scintillators
};

#endif
