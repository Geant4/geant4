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
// $Id: G4SDParticleFilter.hh,v 1.2 2005/11/17 22:53:38 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef G4SDParticleFilter_h
#define G4SDParticleFilter_h 1

class G4Step;
class G4ParticleDefinition;
#include "globals.hh"
#include "G4VSDFilter.hh"

#include <vector>

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//  This class filters steps by partilce definition.
// The particles are given at constructor or add() method.
//
// Created: 2005-11-14  Tsukasa ASO.
// 
///////////////////////////////////////////////////////////////////////////////

class G4SDParticleFilter : public G4VSDFilter 
{

  public: // with description
      G4SDParticleFilter(G4String name);
      G4SDParticleFilter(G4String name,const G4String& particleName);
      G4SDParticleFilter(G4String name,
			 const std::vector<G4String>&  particleNames);
      G4SDParticleFilter(G4String name,
			 const std::vector<G4ParticleDefinition*>&  particleDef);
    // Constructors. Filter name and particle's name.
    //

      virtual ~G4SDParticleFilter();

  public: // with description
      virtual G4bool Accept(const G4Step*) const;

      void add(const G4String& particleName);
      // set method for acceptable particle name.
      //
      void show();

  private:
      std::vector<G4ParticleDefinition*> thePdef;

};

#endif

