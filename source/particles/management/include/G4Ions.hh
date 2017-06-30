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
// $Id: G4Ions.hh 103892 2017-05-03 08:11:00Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 27 June 1998
// ----------------------------------------------------------------
//      Add excitation energy         17 Aug. 1999 H.Kurashige
//      Add isomer level              30 Apr. H.Kurashige


#ifndef G4Ions_h
#define G4Ions_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

// ######################################################################
// ###                          Ions                                 ###
// ######################################################################

class G4Ions : public G4ParticleDefinition
{
 // Class Description
 //  This is the base class for all nuclei including pre-defined 
 //  light nuclei such as deuteron, alpha, and proton (Hydrogen) 
 //  All nuclei/ions created on the fly are objects of this class
 //  Atomic number and atomic mass are vaild only for particles derived
 //  from this class.  This class has Excitation Energy in addition to
 //  the normal particle properties.

 protected:
   G4Ions(){};


 public: //With Description
   G4Ions(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable,  G4bool              shortlived,
       const G4String&     subType ="",
       G4int               anti_encoding =0,
       G4double            excitation = 0.0, 
       G4int               isomer = 0
   );

 public:
   virtual    			~G4Ions();
   G4Ions*    			IonsDefinition();
   G4Ions*    			Ions();

 public:  //With Description
   // Get excitation energy of nucleus
   G4double GetExcitationEnergy() const ; 
  
  // Get Isomer level (=0 for ground state)
  G4int GetIsomerLevel() const; 
   
  // enumerator for floating level base
  enum class G4FloatLevelBase
       { no_Float=0,
         plus_X, plus_Y, plus_Z, plus_U, plus_V, plus_W,
         plus_R, plus_S, plus_T, plus_A, plus_B, plus_C, plus_D, plus_E
       };
  static G4Ions::G4FloatLevelBase FloatLevelBase(char flbChar);
  static G4Ions::G4FloatLevelBase FloatLevelBase(G4int flbIdx);
  static char FloatLevelBaseChar(G4Ions::G4FloatLevelBase flb);

  // set/get methods for floating level base
  G4Ions::G4FloatLevelBase GetFloatLevelBase() const;
  G4int GetFloatLevelBaseIndex() const;
  void SetFloatLevelBase(G4Ions::G4FloatLevelBase flb);
  void SetFloatLevelBase(char flbChar);
  void SetFloatLevelBase(G4int flbIdx);

 private:
  G4double theExcitationEnergy; 
  G4int    theIsomerLevel;
  G4FloatLevelBase floatLevelBase;

};

#define noFloat G4Ions::G4FloatLevelBase::no_Float
#define plusU G4Ions::G4FloatLevelBase::plus_U 
#define plusV G4Ions::G4FloatLevelBase::plus_V 
#define plusW G4Ions::G4FloatLevelBase::plus_W 
#define plusX G4Ions::G4FloatLevelBase::plus_X
#define plusY G4Ions::G4FloatLevelBase::plus_Y 
#define plusZ G4Ions::G4FloatLevelBase::plus_Z 
#define plusR G4Ions::G4FloatLevelBase::plus_R 
#define plusS G4Ions::G4FloatLevelBase::plus_S 
#define plusT G4Ions::G4FloatLevelBase::plus_T 
#define plusA G4Ions::G4FloatLevelBase::plus_A
#define plusB G4Ions::G4FloatLevelBase::plus_B 
#define plusC G4Ions::G4FloatLevelBase::plus_C 
#define plusD G4Ions::G4FloatLevelBase::plus_D 
#define plusE G4Ions::G4FloatLevelBase::plus_E 

inline
 G4Ions* G4Ions::IonsDefinition()
{
  return this;
}

inline
 G4Ions* G4Ions::Ions() 
{
  return this;
}

inline
 G4double G4Ions::GetExcitationEnergy() const 
{
  return theExcitationEnergy;
}

inline
 G4int G4Ions::GetIsomerLevel() const
{
  return theIsomerLevel;
}
    
inline
 G4Ions::G4FloatLevelBase G4Ions::GetFloatLevelBase() const
{
  return floatLevelBase;
}

inline
 G4int G4Ions::GetFloatLevelBaseIndex() const
{
  return static_cast<G4int>(floatLevelBase);
}

inline
 void G4Ions::SetFloatLevelBase(G4Ions::G4FloatLevelBase flb)
{
  floatLevelBase = flb;
}

inline
 void G4Ions::SetFloatLevelBase(char flbChar)
{
  floatLevelBase = FloatLevelBase(flbChar);
}

inline
 void G4Ions::SetFloatLevelBase(G4int flbIdx)
{
  floatLevelBase = FloatLevelBase(flbIdx);
}

inline
 G4Ions::G4FloatLevelBase G4Ions::FloatLevelBase(char flbChar)
{
  G4Ions::G4FloatLevelBase flb = noFloat;
  switch(flbChar)
  {
   case 'x': case 'X':
    flb = plusX;
    break;
   case 'y': case 'Y':
    flb = plusY;
    break;
   case 'z': case 'Z':
    flb = plusZ;
    break;
   case 'u': case 'U':
    flb = plusU;
    break;
   case 'v': case 'V':
    flb = plusV;
    break;
   case 'w': case 'W':
    flb = plusW;
    break;
   case 'r': case 'R':
    flb = plusR;
    break;
   case 's': case 'S':
    flb = plusS;
    break;
   case 't': case 'T':
    flb = plusT;
    break;
   case 'a': case 'A':
    flb = plusA;
    break;
   case 'b': case 'B':
    flb = plusB;
    break;
   case 'c': case 'C':
    flb = plusC;
    break;
   case 'd': case 'D':
    flb = plusD;
    break;
   case 'e': case 'E':
    flb = plusE;
    break;
   case '\0': default:
    break;
  }
  return flb;
}

inline
 G4Ions::G4FloatLevelBase G4Ions::FloatLevelBase(G4int flbIdx)
{
  static G4Ions::G4FloatLevelBase flb[] = 
  { noFloat,
    plusX, plusY, plusZ, plusU, plusV, plusW, 
    plusR, plusS, plusT, plusA, plusB, plusC, plusD, plusE };
  return flb[flbIdx];
}

inline
 char G4Ions::FloatLevelBaseChar(G4Ions::G4FloatLevelBase flb)
{
  static char flbChar[] = {'\0','X','Y','Z','U','V','W',
                                'R','S','T','A','B','C','D','E'};
  return flbChar[static_cast<G4int>(flb)];
}

#endif








