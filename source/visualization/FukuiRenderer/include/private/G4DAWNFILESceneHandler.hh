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
//
//
// Satoshi TANAKA

#ifndef G4DAWNFILE_SCENE_HANDLER_HH
#define G4DAWNFILE_SCENE_HANDLER_HH

#include "globals.hh"

#include "G4VSceneHandler.hh"

#include "G4FRofstream.hh"
#include "G4FRConst.hh"

class G4VisAttributes;
class G4DAWNFILE;

//-----
class G4DAWNFILESceneHandler : public G4VSceneHandler
{
  friend class G4DAWNFILEViewer;

 public:
  //----- constructor and destructor
  G4DAWNFILESceneHandler(G4DAWNFILE& system, const G4String& name = "");
  virtual ~G4DAWNFILESceneHandler();

  //----- overriding base class methods
  using G4VSceneHandler::AddPrimitive;
  void AddPrimitive(const G4Polyline& line);
  void AddPrimitive(const G4Polyhedron& p);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);

  void BeginPrimitives(const G4Transform3D& objectTransformation);
  void EndPrimitives();

  using G4VSceneHandler::AddSolid;
  void AddSolid(const G4Box& box);
  void AddSolid(const G4Cons& cons);
  void AddSolid(const G4Tubs& tubs);
  void AddSolid(const G4Trd& trd);
  void AddSolid(const G4Trap& trap);
  void AddSolid(const G4Sphere& sphere);
  void AddSolid(const G4Para& para);
  void AddSolid(const G4Torus& torus);
  void AddSolid(const G4VSolid& solid);

  void ClearTransientStore();  // Used for triggering detector re-drawing.

  //----- public methods inherent to this class
  void FRBeginModeling();
  void FREndModeling();
  G4bool FRIsInModeling() { return FRflag_in_modeling; }

  G4bool IsSavingG4Prim(void) { return flag_saving_g4_prim; }
  void BeginSavingG4Prim(void);
  void EndSavingG4Prim(void);
  void SetG4PrimFileName();

  G4DAWNFILE& GetSystem() { return fSystem; }
  void SendBoundingBox(void);
  const char* GetG4PrimFileName() { return fG4PrimFileName; }

 private:
  //----- Utilities etc (common to DAWN and DAWNFILE drivers )
  G4bool SendVisAttributes(const G4VisAttributes* pAV);
  G4bool IsVisible(void);
  void SendTransformedCoordinates(void);
  void SendPhysVolName(void);
  void SendNdiv(void);

  //----- public methods common to DAWN and DAWNFILE drivers
 public:
  void SendStr(const char* char_string);
  void SendStrInt(const char* char_string, G4int ival);
  void SendStrInt3(const char* char_string, G4int ival1, G4int ival2,
                   G4int ival3);
  void SendStrInt4(const char* char_string, G4int ival1, G4int ival2,
                   G4int ival3, G4int ival4);
  void SendStrDouble(const char* char_string, G4double dval);
  void SendStrDouble2(const char* char_string, G4double dval1, G4double dval2);
  void SendStrDouble3(const char* char_string, G4double dval1, G4double dval2,
                      G4double dval3);

  void SendStrDouble4(const char* char_string, G4double dval1, G4double dval2,
                      G4double dval3, G4double dval4);

  void SendStrDouble5(const char* char_string, G4double dval1, G4double dval2,
                      G4double dval3, G4double dval4, G4double dval5);

  void SendStrDouble6(const char* char_string, G4double dval1, G4double dval2,
                      G4double dval3, G4double dval4, G4double dval5,
                      G4double dval6);

  void SendStrDouble7(const char* char_string, G4double dval1, G4double dval2,
                      G4double dval3, G4double dval4, G4double dval5,
                      G4double dval6, G4double dval7);

  void SendStrDouble11(const char* char_string, G4double dval1, G4double dval2,
                       G4double dval3, G4double dval4, G4double dval5,
                       G4double dval6, G4double dval7, G4double dval8,
                       G4double dval9, G4double dval10, G4double dval11);

  void SendIntDouble3(G4int ival, G4double dval1, G4double dval2,
                      G4double dval3);
  void SendInt3Str(G4int ival1, G4int ival2, G4int ival3,
                   const char* char_string);
  void SendInt4Str(G4int ival1, G4int ival2, G4int ival3, G4int ival4,
                   const char* char_string);

  void SendStrDouble3Str(const char* char_string1, G4double dval1,
                         G4double dval2, G4double dval3,
                         const char* char_string2);

  void SendStrDouble6Str(const char* char_string1, G4double dval1,
                         G4double dval2, G4double dval3, G4double dval4,
                         G4double dval5, G4double dval6,
                         const char* char_string2);

  void SendInt(G4int val);
  void SendDouble(G4double val);

 private:
  G4DAWNFILE& fSystem;  // Graphics system for this scene.
  static G4int fSceneIdCount;

  G4FRofstream fPrimDest;  // defined here
  G4bool FRflag_in_modeling;
  // true:  FR_BEGIN_MODELING has sent to DAWN, and
  //        FR_END_MODELING   has not sent yet.
  // false:  otherwise
  //
  // FRflag_in_modeling is set to "true"
  // in FRBeginModeling(), and to "false"
  // in FREndModeling().

  G4bool flag_saving_g4_prim;

  const int COMMAND_BUF_SIZE;

  char fG4PrimDestDir[256];
  char fG4PrimFileName[256];
  G4int fMaxFileNum;

  G4int fPrec, fPrec2;
};

#endif
