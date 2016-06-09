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
// $Id: MedLinacMLCDecorator.hh,v 1.2 2004/11/25 16:28:39 mpiergen Exp $
//
// Code developed by: M. Piergentili
//
//


#ifndef MedLinacMLCDecorator_h
#define MedLinacMLCDecorator_h 1

#include "globals.hh"
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

  void SetLeafName (G4String);
  void SetPos_y (G4double);

private:
  void ConstructMultiLeafCollimator(G4VPhysicalVolume*,G4VPhysicalVolume*);

   G4LogicalVolume* leafLog;
   G4VPhysicalVolume* leafAPhys;
   G4VPhysicalVolume* leafBPhys;

  G4VisAttributes* simpleTungstenSVisAtt;

  public:
  
  void PrintParametersMLC(); 

  G4String GetLeafName()  {return leaf_name;}; 
  G4double GetPos_y()  {return pos;}; 
 


  private:
  G4String  leaf_name;
  G4double  pos;

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

