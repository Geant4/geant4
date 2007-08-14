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
// $Id: G4SDProcessFilter.hh,v 1.1 2007-08-14 21:23:51 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4SDProcessFilter_h
#define G4SDProcessFilter_h 1

class G4Step;
class G4VProcess;
class G4ParticleDefinition;
#include "globals.hh"
#include "G4VSDFilter.hh"

#include <vector>

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector. 
//  This class filters steps by process of partilce.
// The acceptable process are given at constructor or add() method.
//
// Created: 2007-03-22  Tsukasa ASO.
// Modified: 2007-08-14 T.Aso Process should be given as an G4VProcess Object
//                      rather than processName because of dependency of 
//                      category. 
// 
///////////////////////////////////////////////////////////////////////////////

class G4SDProcessFilter : public G4VSDFilter 
{

  public: // with description
      G4SDProcessFilter(G4String name);
      G4SDProcessFilter(G4String name,
			const G4VProcess* process, const G4String& particleName );
    // Constructors. Filter name and process's name.
    //

      virtual ~G4SDProcessFilter();
  // Destrcutor

  public:
     void add(const G4VProcess* process, const G4String& particleName); 
  
  public: // with description
      virtual G4bool Accept(const G4Step*) const;

      void show();

  private:
     typedef std::vector<const G4VProcess*>  G4SDProcessCollection; 
     G4SDProcessCollection  theProcessDef;
     typedef std::vector<G4ParticleDefinition*> G4SDParticleCollection;
     G4SDParticleCollection theParticleDef;
};

#endif

