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















