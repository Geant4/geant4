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
#ifndef G4QGSModel_h
#define G4QGSModel_h 1

// Class Description
// Model for hadron (p,n,pi,K) nuclear reactions in geant4. IT implements
// A. Kaydalov's quark gluon string model.
// To be used in your physics list, in case you need this kind of physics.
// Class Description - End

#include <cmath>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ExcitedStringVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4VPartonStringModel.hh"
#include "G4QGSParticipants.hh"
#include "G4DiffractiveStringBuilder.hh"
#include "G4SoftStringBuilder.hh"
#include "G4PartonPair.hh"

//*********************************************************************************************** 


//*****************************************************************************************

template<class ParticipantType>
class G4QGSModel : public G4VPartonStringModel
    {
// Constructors   
public:
    G4QGSModel();
    virtual ~G4QGSModel();
    G4QGSModel(const G4QGSModel &right);
    G4QGSModel& operator=(const G4QGSModel &right);

// Method
public:
    virtual G4V3DNucleus* GetWoundedNucleus() const;
    virtual G4V3DNucleus* GetProjectileNucleus() const;  // Uzhi Nov. 2012
    virtual void Init(const G4Nucleus& Nucleus, const G4DynamicParticle& Projectile);
    virtual G4ExcitedStringVector * GetStrings();
    virtual void ModelDescription(std::ostream& outFile) const;

private:
   ParticipantType theParticipants;
   G4DiffractiveStringBuilder theDiffractiveStringBuilder;
   G4SoftStringBuilder theSoftStringBuilder;

   };

//-------------------------------------------------------------------------------------------


//*****************************************************************************************
    
#include "G4QGSModel.icc"

#endif


