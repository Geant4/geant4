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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2VParticleReporter.hh,v 1.3 2001-07-11 10:02:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2VParticleReporter_h
#define tst2VParticleReporter_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "tst2ParticleContainer.hh"

class tst2VParticleReporter
{
 public:
  //constructors
    tst2VParticleReporter();

  //destructor
    virtual ~tst2VParticleReporter();

 public:
    // equality operators
    G4int operator==(const tst2VParticleReporter &right) const 
    {   return (this == &right);    }

    G4int operator!=(const tst2VParticleReporter &right) const 
    {   return (this != &right);    }

 public:
	virtual void Print(const tst2ParticleContainer& container, 
                       const G4String& option) = 0;

 protected:
    G4ParticleDefinition* GetParticle(G4int index) const;
    G4int                 GetEncoding(G4int index) const;
    G4ParticleDefinition* GetAntiParticle(G4int index) const;
    G4int                 entries() const;
    const tst2ParticleContainer*  pList;
  
};

inline     
 G4int tst2VParticleReporter::entries() const
{
  return pList->entries();
}

inline     
 G4ParticleDefinition*  tst2VParticleReporter::GetParticle(G4int index) const
{
  return pList->GetParticle(index);
}

inline     
 G4ParticleDefinition*  tst2VParticleReporter::GetAntiParticle(G4int index) const
{
  return pList->GetAntiParticle(index);
}

inline     
 G4int tst2VParticleReporter::GetEncoding(G4int index) const
{
  return pList->GetEncoding(index);
}
 


#endif
