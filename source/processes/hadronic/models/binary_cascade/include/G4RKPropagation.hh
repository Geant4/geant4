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

#ifndef G4RKPropagation_h
#define G4RKPropagation_h 1

#include "G4VFieldPropagation.hh"
#include "G4VNuclearField.hh"
#include "G4V3DNucleus.hh"
#include "G4KM_DummyField.hh"
#include "G4Mag_EqRhs.hh"
#include <map>


class G4RKPropagation: public G4VFieldPropagation
{

public:
  G4RKPropagation();
  virtual ~G4RKPropagation();

private:
  G4RKPropagation(const  G4RKPropagation &right);
  const G4RKPropagation & operator=(const G4RKPropagation & right);
  G4int operator==(const G4RKPropagation & right) const;
  G4int operator!=(const G4RKPropagation & right) const;

public:

  virtual void Init(G4V3DNucleus * nucleus);
  virtual void Transport(G4KineticTrackVector &theActive,
			 const G4KineticTrackVector &theSpectators,
			 G4double theTimeStep);
  G4bool GetSphereIntersectionTimes(const G4KineticTrack * track,
				    G4double & t1, G4double & t2);
  G4ThreeVector GetMomentumTransfer() const;
private:
  G4double theOuterRadius;
  G4V3DNucleus * theNucleus;
  std::map <G4int, G4VNuclearField *, std::less<G4int> > * theFieldMap;
  std::map <G4int, G4Mag_EqRhs *, std::less<G4int> > * theEquationMap;
  G4KM_DummyField * theField;

  G4ThreeVector theMomentumTranfer;

  G4bool GetSphereIntersectionTimes(const G4double radius,
				    const G4ThreeVector & currentPos,
				    const G4LorentzVector & momentum,
				    G4double & t1, G4double & t2);
// implementation

  G4bool FieldTransport(G4KineticTrack * track, const G4double timestep);
  G4bool FreeTransport(G4KineticTrack * track, const G4double timestep);

  void delete_FieldsAndMap(
	std::map <G4int, G4VNuclearField *, std::less<G4int> > * aMap);
  void delete_EquationsAndMap(
	std::map <G4int, G4Mag_EqRhs *, std::less<G4int> > * aMap);
  
  public:
  inline G4double GetBarrier(G4int encoding) 
  {     
	std::map <G4int, G4VNuclearField *, std::less<G4int> >::iterator iter;
	iter = theFieldMap->find(encoding);
	if(iter == theFieldMap->end()) return 0;
	return (*theFieldMap)[encoding]->GetBarrier();
  }
  
  inline G4double GetField(G4int encoding,G4ThreeVector pos)
  {
	std::map <G4int, G4VNuclearField *, std::less<G4int> >::iterator iter;
	iter = theFieldMap->find(encoding);
	if(iter == theFieldMap->end()) return 0;
	return (*theFieldMap)[encoding]->GetField(pos);
  }

};

#endif















