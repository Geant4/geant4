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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4XnpElasticLowE
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4XNPELASTICLOWE_HH
#define G4XNPELASTICLOWE_HH

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"
#include "G4PhysicsVector.hh"

class G4KineticTrack;

class G4XnpElasticLowE : public G4VCrossSectionSource
{

public:

  G4XnpElasticLowE();

  virtual ~G4XnpElasticLowE();

  G4bool operator==(const G4XnpElasticLowE &right) const;
  G4bool operator!=(const G4XnpElasticLowE &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4bool IsValid(G4double e) const;

  virtual void Print() const;

  virtual G4String Name() const;

  virtual G4double HighLimit() const { return _highLimit; }


protected:


private:  

  G4XnpElasticLowE(const G4XnpElasticLowE &right);
  const G4XnpElasticLowE& operator=(const G4XnpElasticLowE &right);
  
  static const G4double _lowLimit;
  static const G4double _highLimit;
  static const G4double _sigmaTable[101];
  static const G4int _tableSize;
  static const G4double _eMinTable;
  static const G4double _eStepLog;

  G4PhysicsVector* _sigma;
  G4double _eMin;
  G4double _eMax;

};

#endif











































