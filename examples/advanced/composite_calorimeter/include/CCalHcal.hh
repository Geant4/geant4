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

class CCalHcal: public CCalDetector {
public:
  //Constructor and Destructor
  CCalHcal(const G4String &name);
  virtual ~CCalHcal();

  //Get Methods
  G4String getGenMat()                    const {return genMaterial;}
  double   getDy_2Cal()                   const {return dy_2Cal;}
  double   getDx_2Cal()                   const {return dx_2Cal;}
  double   getXposCal()                   const {return xposCal;}
  G4String getBoxMat()                    const {return boxMaterial;}
  int      getNBox()                      const {return nBox;}
  double   getDy_2Box()                   const {return dy_2Box;}
  double   getDx_2Box()                   const {return dx_2Box;}
  double   getWallThickBox()              const {return wallThickBox;}
  double   getXposBox(unsigned int i)     const {return xposBox[i];}
  int      getNLayerScnt()                const {return nLayerScnt;}
  int      getTypeScnt(unsigned int i)    const {return typeLayerScnt[i];}
  int      getMotherScnt(unsigned int i)  const {return mothLayerScnt[i];}
  double   getXposScnt(unsigned int i)    const {return xposLayerScnt[i];}
  int      getNLayerAbs()                 const {return nLayerAbs;}
  int      getTypeAbs(unsigned int i)     const {return typeLayerAbs[i];}
  int      getMotherAbs(unsigned int i)   const {return mothLayerAbs[i];}
  double   getXposAbs(unsigned int i)     const {return xposLayerAbs[i];}
  G4String getAbsMat()                    const {return absMaterial;}
  int      getNAbsorber()                 const {return nAbsorber;}
  double   getDy_2Abs(     )              const {return dy_2Absorber;}
  double   getDx_2Abs(unsigned int i)     const {return dx_2Absorber[i];}
  G4String getScntMat()                   const {return scntMaterial;}
  G4String getWrapMat()                   const {return wrapMaterial;}
  G4String getPlasMat()                   const {return plasMaterial;}
  int      getNScintillator()             const {return nScintillator;}
  double   getDy_2ScntLay(unsigned int i) const {return dy_2ScntLayer[i];}
  double   getDx_2ScntLay(unsigned int i) const {return dx_2ScntLayer[i];}
  double   getDx_2Wrap(unsigned int i)    const {return dx_2Wrapper[i];}
  double   getDx_2FrontP(unsigned int i)  const {return dx_2FrontPlastic[i];}
  double   getDx_2BackP(unsigned int i)   const {return dx_2BackPlastic[i];}
  double   getDx_2Scnt(unsigned int i)    const {return dx_2Scintillator[i];}

protected:
  virtual int readFile();
  virtual void constructDaughters();

private:
  G4String genMaterial;            //General material
  double   dy_2Cal;                //Half width     of the Hcal 
  double   dx_2Cal;                //Half thickness of the Hcal
  double   xposCal;                //Position in mother

  G4String boxMaterial;            //Material of boxes
  int      nBox;                   //Number of boxes
  double   dy_2Box;                //Half width     of the Boxes
  double   dx_2Box;                //Half thickness of the Boxes
  double   wallThickBox;           //Wall thickness of the boxes
  double*  xposBox;                //Position in mother

  int      nLayerScnt;             //Number of scintillator layers
  int*     typeLayerScnt;          //Layer type
  int*     mothLayerScnt;          //Mother type
  double*  xposLayerScnt;          //Position in mother

  int      nLayerAbs;              //Number of absorber     layers
  int*     typeLayerAbs;           //Layer type
  int*     mothLayerAbs;           //Mother type
  double*  xposLayerAbs;           //Position in mother

  G4String absMaterial;            //Material of absorbers
  int      nAbsorber;              //Number of absorber types
  double   dy_2Absorber;           //Half width     of the absorbers
  double*  dx_2Absorber;           //Half thickness of the absorbers

  G4String scntMaterial;           //Material of Scintillator
  G4String wrapMaterial;           //Material of Wrapper
  G4String plasMaterial;           //Material of plastic cover
  int      nScintillator;          //Number of scintillator types
  double*  dy_2ScntLayer;          //Half width     of scintillator layers
  double*  dx_2ScntLayer;          //Half thickness of scintillator layers
  double*  dx_2Wrapper;            //Half thickness of wrappers
  double*  dx_2FrontPlastic;       //Half thickness of front plastic
  double*  dx_2BackPlastic;        //Half thickness of back  plastic
  double*  dx_2Scintillator;       //Half thickness of scintillators
};

#endif
