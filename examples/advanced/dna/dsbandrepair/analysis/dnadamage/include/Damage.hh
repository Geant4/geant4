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
/// \file Damage.hh
/// \brief Definition of the Damage class

#ifndef DAMAGE_HH
#define DAMAGE_HH

struct Position
{
  Position(double v1,double v2,double v3): x(v1), y(v2), z(v3) {}
  double x=0;
  double y=0;
  double z=0;
  void setX(double v) {x=v;}
  void setY(double v) {y=v;}
  void setZ(double v) {z=v;}
};

/// \brief defines a DNA damage
class Damage
{
public:

/** DNA Damage type
*   Defines the molecule type of damaged DNA following SDD formalism 
*/
  enum DamageType{
    fOther = -1,
    fBackbone = 0,
    fBase = 1
  };

/** DNA Damage cause
*   Defines the cause of DNA damage following SDD formalism 
*/
  enum DamageCause{
    fUnknown = -1,
    fDirect = 0,
    fIndirect = 1
  };

/** Damaged DNA
*   Defines the damaged DNA structure following SDD formalism 
*/
  enum DamageChromatin{
    fUnspecified = 0,
    fHeterochromatin = 1,
    fEuchromatin = 2,
    fFreeDNA = 3,
    fOtherDNA = 4
  };
  
  /// \brief constructor
  Damage(DamageType,unsigned int,unsigned int,unsigned int,unsigned long int,Position,DamageCause,DamageChromatin);
  /// \brief destructor
  ~Damage() = default;

  // Getters and setters

  void SetDamageType(DamageType pVal){fType=pVal;};
  DamageType GetDamageType(){return fType;};

  void SetChromo(unsigned int pVal){fChromo=pVal;};
  unsigned int GetChromo() const{return fChromo;};

  void SetEvt(unsigned int pVal){fEvt=pVal;};
  unsigned int GetEvt() const{return fEvt;};

  void SetStrand(unsigned int pVal){fStrand=pVal;};
  unsigned int GetStrand() const{return fStrand;};

  void SetCopyNb(unsigned long int pVal){fCopyNb=pVal;};
  unsigned long int GetCopyNb() const{return fCopyNb;};

  void SetCause(DamageCause pVal){fCause=pVal;};
  DamageCause GetCause() const{return fCause;};

  void SetDamageChromatin(DamageChromatin pVal){fChromatin=pVal;};
  DamageChromatin GetDamageChromatin(){return fChromatin;};

  bool operator != (const Damage& ) const;
  bool operator == (const Damage& ) const;

private:
  
  DamageType        fType;    // SB or BD?
  unsigned int      fChromo;  // chromosome ID
  unsigned int      fEvt;     // event number
  unsigned int      fStrand;  // Strand
  unsigned long int fCopyNb;  // Copy number
  Position             fPosition;// Position
  DamageCause       fCause;  // Direct or indirect damage?
  DamageChromatin   fChromatin;   // hetero or euchromatin?
};

#endif // DAMAGE_HH