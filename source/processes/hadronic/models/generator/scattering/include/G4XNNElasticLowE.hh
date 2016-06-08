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

#ifndef G4XNNElasticLowE_h
#define G4XNNElasticLowE_h

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4CrossSectionVector.hh"
#include "G4PhysicsVector.hh"
#include "g4std/map"

class G4KineticTrack;

class G4XNNElasticLowE : public G4VCrossSectionSource
{

public:

  G4XNNElasticLowE();

  virtual ~G4XNNElasticLowE();

  G4bool operator==(const G4XNNElasticLowE &right) const;
  G4bool operator!=(const G4XNNElasticLowE &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;
 
  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4bool IsValid(G4double e) const;

  virtual G4String Name() const;

  virtual void Print() const;

  virtual G4double HighLimit() const { return _highLimit; }


protected:
 

private:  

  G4XNNElasticLowE(const G4XNNElasticLowE &right);
  const G4XNNElasticLowE& operator=(const G4XNNElasticLowE &right);
  
  static const G4double _lowLimit;
  static const G4double _highLimit;
  static const G4double ppTable[101];
  static const G4double npTable[101];
  static const G4int tableSize;
  static const G4double _eMinTable;
  static const G4double _eStepLog;

  G4std::map <G4String, G4PhysicsVector*, G4std::less<G4String> > xMap;

  G4double _eMin;
  G4double _eMax;

};

#endif











































