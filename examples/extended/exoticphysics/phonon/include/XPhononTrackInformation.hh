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
/// \file exoticphysics/phonon/include/XPhononTrackInformation.hh
/// \brief Definition of the XPhononTrackInformation class
//
// $Id$
//
#ifndef XPhononTrackInformation_h
#define XPhononTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class XPhononTrackInformation : public G4VUserTrackInformation
{

public:
  XPhononTrackInformation();
  XPhononTrackInformation(const G4Track* aTrack);
  XPhononTrackInformation(const XPhononTrackInformation* aTrackInfo);
  XPhononTrackInformation(G4ThreeVector kNew);
  virtual ~XPhononTrackInformation();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator == (const XPhononTrackInformation& right) const {return (this==&right);}

  XPhononTrackInformation& operator = (const XPhononTrackInformation& right);

  //void SetSourceTrackInformation(const G4Track* aTrack);
  void Print() const;

private:
  G4ThreeVector fK;

public:
  inline G4ThreeVector GetK() const {return fK;}
  inline void SetK(G4ThreeVector kNew){fK=G4ThreeVector(kNew);}

};

extern G4Allocator<XPhononTrackInformation> aTrackInformationAllocator;

inline void* XPhononTrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;

}

inline void XPhononTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((XPhononTrackInformation*)aTrackInfo);}

#endif
