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
// $Id: MedLinacMLCDecorator.hh,v 1.1 2004-11-24 16:53:29 mpiergen Exp $
//
// Code developed by: M. Piergentili
//
//


#ifndef MedLinacMLCDecorator_h
#define MedLinacMLCDecorator_h 1

#include "MedLinacDecorator.hh"
#include "MedLinacVGeometryComponent.hh"
#include "MedLinacMLCMessenger.hh"

class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class G4VisAttributes;
class MedLinacVGeometryComponent;
class MedLinacDecorator;
class MedLinacMLCMessenger;

class MedLinacMLCDecorator: public MedLinacDecorator
{
public:
  MedLinacMLCDecorator(MedLinacVGeometryComponent*);
  ~MedLinacMLCDecorator();
  void ConstructComponent(G4VPhysicalVolume*,G4VPhysicalVolume*);
  void DestroyComponent(); 

  void SetA1Pos_y (G4double);
  void SetA2Pos_y (G4double);
  void SetA3Pos_y (G4double);
  void SetA4Pos_y (G4double);
  void SetA5Pos_y (G4double);
  void SetA6Pos_y (G4double);
  void SetA7Pos_y (G4double);
  void SetA8Pos_y (G4double);
  void SetA9Pos_y (G4double);
  void SetA10Pos_y (G4double);
  void SetA11Pos_y (G4double);
  void SetA12Pos_y (G4double);
  void SetA13Pos_y (G4double);
  void SetA14Pos_y (G4double);
  void SetA15Pos_y (G4double);
  void SetA16Pos_y (G4double);
  void SetA17Pos_y (G4double);
  void SetA18Pos_y (G4double);
  void SetA19Pos_y (G4double);
  void SetA20Pos_y (G4double);
  void SetA21Pos_y (G4double);
  void SetA22Pos_y (G4double);
  void SetA23Pos_y (G4double);
  void SetA24Pos_y (G4double);
  void SetA25Pos_y (G4double);
  void SetA26Pos_y (G4double);
  void SetA27Pos_y (G4double);
  void SetA28Pos_y (G4double);
  void SetA29Pos_y (G4double);
  void SetA30Pos_y (G4double);
  void SetA31Pos_y (G4double);
  void SetA32Pos_y (G4double);
  void SetA33Pos_y (G4double);
  void SetA34Pos_y (G4double);
  void SetA35Pos_y (G4double);
  void SetA36Pos_y (G4double);
  void SetA37Pos_y (G4double);
  void SetA38Pos_y (G4double);
  void SetA39Pos_y (G4double);
  void SetA40Pos_y (G4double);


  void SetB1Pos_y (G4double);
  void SetB2Pos_y (G4double);
  void SetB3Pos_y (G4double);
  void SetB4Pos_y (G4double);
  void SetB5Pos_y (G4double);
  void SetB6Pos_y (G4double);
  void SetB7Pos_y (G4double);
  void SetB8Pos_y (G4double);
  void SetB9Pos_y (G4double);
  void SetB10Pos_y (G4double);
  void SetB11Pos_y (G4double);
  void SetB12Pos_y (G4double);
  void SetB13Pos_y (G4double);
  void SetB14Pos_y (G4double);
  void SetB15Pos_y (G4double);
  void SetB16Pos_y (G4double);
  void SetB17Pos_y (G4double);
  void SetB18Pos_y (G4double);
  void SetB19Pos_y (G4double);
  void SetB20Pos_y (G4double);
  void SetB21Pos_y (G4double);
  void SetB22Pos_y (G4double);
  void SetB23Pos_y (G4double);
  void SetB24Pos_y (G4double);
  void SetB25Pos_y (G4double);
  void SetB26Pos_y (G4double);
  void SetB27Pos_y (G4double);
  void SetB28Pos_y (G4double);
  void SetB29Pos_y (G4double);
  void SetB30Pos_y (G4double);
  void SetB31Pos_y (G4double);
  void SetB32Pos_y (G4double);
  void SetB33Pos_y (G4double);
  void SetB34Pos_y (G4double);
  void SetB35Pos_y (G4double);
  void SetB36Pos_y (G4double);
  void SetB37Pos_y (G4double);
  void SetB38Pos_y (G4double);
  void SetB39Pos_y (G4double);
  void SetB40Pos_y (G4double);

private:
  void ConstructMultiLeafCollimator(G4VPhysicalVolume*,G4VPhysicalVolume*);

   G4LogicalVolume* leafLog;
   G4VPhysicalVolume* leafAPhys;
   G4VPhysicalVolume* leafBPhys;

  G4VisAttributes* simpleTungstenSVisAtt;

  public:
  
   void PrintParametersMLC(); 

  G4double GetA1Pos_y()  {return a1y;}; 
  G4double GetA2Pos_y()  {return a2y;}; 
  G4double GetA3Pos_y()  {return a3y;}; 
  G4double GetA4Pos_y()  {return a4y;}; 
  G4double GetA5Pos_y()  {return a5y;}; 
  G4double GetA6Pos_y()  {return a6y;}; 
  G4double GetA7Pos_y()  {return a7y;}; 
  G4double GetA8Pos_y()  {return a8y;}; 
  G4double GetA9Pos_y()  {return a9y;}; 
  G4double GetA10Pos_y()  {return a10y;}; 
  G4double GetA11Pos_y()  {return a11y;}; 
  G4double GetA12Pos_y()  {return a12y;}; 
  G4double GetA13Pos_y()  {return a13y;}; 
  G4double GetA14Pos_y()  {return a14y;}; 
  G4double GetA15Pos_y()  {return a15y;}; 
  G4double GetA16Pos_y()  {return a16y;}; 
  G4double GetA17Pos_y()  {return a17y;}; 
  G4double GetA18Pos_y()  {return a18y;}; 
  G4double GetA19Pos_y()  {return a19y;}; 
  G4double GetA20Pos_y()  {return a20y;}; 
  G4double GetA21Pos_y()  {return a21y;}; 
  G4double GetA22Pos_y()  {return a22y;}; 
  G4double GetA23Pos_y()  {return a23y;}; 
  G4double GetA24Pos_y()  {return a24y;}; 
  G4double GetA25Pos_y()  {return a25y;}; 
  G4double GetA26Pos_y()  {return a26y;}; 
  G4double GetA27Pos_y()  {return a27y;}; 
  G4double GetA28Pos_y()  {return a28y;}; 
  G4double GetA29Pos_y()  {return a29y;}; 
  G4double GetA30Pos_y()  {return a30y;}; 
  G4double GetA31Pos_y()  {return a31y;}; 
  G4double GetA32Pos_y()  {return a32y;}; 
  G4double GetA33Pos_y()  {return a33y;}; 
  G4double GetA34Pos_y()  {return a34y;}; 
  G4double GetA35Pos_y()  {return a35y;}; 
  G4double GetA36Pos_y()  {return a36y;}; 
  G4double GetA37Pos_y()  {return a37y;}; 
  G4double GetA38Pos_y()  {return a38y;}; 
  G4double GetA39Pos_y()  {return a39y;}; 
  G4double GetA40Pos_y()  {return a40y;}; 

  G4double GetB1Pos_y()  {return b1y;}; 
  G4double GetB2Pos_y()  {return b2y;}; 
  G4double GetB3Pos_y()  {return b3y;}; 
  G4double GetB4Pos_y()  {return b4y;}; 
  G4double GetB5Pos_y()  {return b5y;}; 
  G4double GetB6Pos_y()  {return b6y;}; 
  G4double GetB7Pos_y()  {return b7y;}; 
  G4double GetB8Pos_y()  {return b8y;}; 
  G4double GetB9Pos_y()  {return b9y;}; 
  G4double GetB10Pos_y()  {return b10y;}; 
  G4double GetB11Pos_y()  {return b11y;}; 
  G4double GetB12Pos_y()  {return b12y;}; 
  G4double GetB13Pos_y()  {return b13y;}; 
  G4double GetB14Pos_y()  {return b14y;}; 
  G4double GetB15Pos_y()  {return b15y;}; 
  G4double GetB16Pos_y()  {return b16y;}; 
  G4double GetB17Pos_y()  {return b17y;}; 
  G4double GetB18Pos_y()  {return b18y;}; 
  G4double GetB19Pos_y()  {return b19y;}; 
  G4double GetB20Pos_y()  {return b20y;}; 
  G4double GetB21Pos_y()  {return b21y;}; 
  G4double GetB22Pos_y()  {return b22y;}; 
  G4double GetB23Pos_y()  {return b23y;}; 
  G4double GetB24Pos_y()  {return b24y;}; 
  G4double GetB25Pos_y()  {return b25y;}; 
  G4double GetB26Pos_y()  {return b26y;}; 
  G4double GetB27Pos_y()  {return b27y;}; 
  G4double GetB28Pos_y()  {return b28y;}; 
  G4double GetB29Pos_y()  {return b29y;}; 
  G4double GetB30Pos_y()  {return b30y;}; 
  G4double GetB31Pos_y()  {return b31y;}; 
  G4double GetB32Pos_y()  {return b32y;}; 
  G4double GetB33Pos_y()  {return b33y;}; 
  G4double GetB34Pos_y()  {return b34y;}; 
  G4double GetB35Pos_y()  {return b35y;}; 
  G4double GetB36Pos_y()  {return b36y;}; 
  G4double GetB37Pos_y()  {return b37y;}; 
  G4double GetB38Pos_y()  {return b38y;}; 
  G4double GetB39Pos_y()  {return b39y;}; 
  G4double GetB40Pos_y()  {return b40y;}; 


  private:
  G4double  a1y;
  G4double  a2y;
  G4double  a3y;
  G4double  a4y;
  G4double  a5y;
  G4double  a6y;
  G4double  a7y;
  G4double  a8y;
  G4double  a9y;
  G4double  a10y;
  G4double  a11y;
  G4double  a12y;
  G4double  a13y;
  G4double  a14y;
  G4double  a15y;
  G4double  a16y;
  G4double  a17y;
  G4double  a18y;
  G4double  a19y;
  G4double  a20y;
  G4double  a21y;
  G4double  a22y;
  G4double  a23y;
  G4double  a24y;
  G4double  a25y;
  G4double  a26y;
  G4double  a27y;
  G4double  a28y;
  G4double  a29y;
  G4double  a30y;
  G4double  a31y;
  G4double  a32y;
  G4double  a33y;
  G4double  a34y;
  G4double  a35y;
  G4double  a36y;
  G4double  a37y;
  G4double  a38y;
  G4double  a39y;
  G4double  a40y;

  G4double  b1y;
  G4double  b2y;
  G4double  b3y;
  G4double  b4y;
  G4double  b5y;
  G4double  b6y;
  G4double  b7y;
  G4double  b8y;
  G4double  b9y;
  G4double  b10y;
  G4double  b11y;
  G4double  b12y;
  G4double  b13y;
  G4double  b14y;
  G4double  b15y;
  G4double  b16y;
  G4double  b17y;
  G4double  b18y;
  G4double  b19y;
  G4double  b20y;
  G4double  b21y;
  G4double  b22y;
  G4double  b23y;
  G4double  b24y;
  G4double  b25y;
  G4double  b26y;
  G4double  b27y;
  G4double  b28y;
  G4double  b29y;
  G4double  b30y;
  G4double  b31y;
  G4double  b32y;
  G4double  b33y;
  G4double  b34y;
  G4double  b35y;
  G4double  b36y;
  G4double  b37y;
  G4double  b38y;
  G4double  b39y;
  G4double  b40y;

  MedLinacMLCMessenger* MLCMessenger;
};
#endif

