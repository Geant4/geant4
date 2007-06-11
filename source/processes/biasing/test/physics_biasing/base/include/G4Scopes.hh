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
// $Id: G4Scopes.hh,v 1.2 2007-06-11 19:25:47 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, May 2007. Creation - scope definitions.
//
#ifndef G4SCOPES_HH
#define G4SCOPES_HH
   
#include "G4FunctorIdentifier.hh"
#include "G4ObserverCollectionT.hh"
#include "G4TypeList.hh"
class G4Track;
class G4Step;
    
namespace G4Scopes {
    
  namespace {
  
    G4bool defaultA()
    {
      return true;
    }
    
    G4bool defaultB(const G4Track&)
    {
      return true;
    }
      
    G4bool defaultC(const G4Track&, const G4Step&)
    {
      return true;
    }
  }
      
  template <typename Args> struct Scope
  {
    typedef G4Functor<G4bool, G4FunctorIdentifier, Args> ActivatorDef;
    typedef G4ObserverCollectionT<Args> ObserverCollection;
  };
  struct Global {
    typedef G4NullType Arg;
    typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;  
    typedef G4bool (*ActivatorFunc)();
    typedef G4ObserverCollectionT<Arg> ObserverCollection;
  };

  struct Tracking {
    struct StartTracking {
      typedef G4TypeList_1(G4Track*) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(G4Track*);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };

    struct EndTracking {
      typedef G4TypeList_1(G4Track) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(const G4Track&);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };
  };
     
  struct Stepping {
    struct StartStep {
      typedef G4TypeList_2(G4Track, G4Step) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(const G4Track&, const G4Step&);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };
      
    struct EndStep {
      typedef G4TypeList_2(G4Track, G4Step) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(const G4Track&, const G4Step&);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };
  };
  struct GeometryBoundary {
    struct StartBoundary {
      typedef G4TypeList_2(G4Track, G4Step) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(const G4Track&, const G4Step&);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };
      
    struct EndBoundary {
      typedef G4TypeList_2(G4Track, G4Step) Arg;
      typedef G4Functor<G4bool, G4FunctorIdentifier, Arg> ActivatorDef;
      typedef G4bool (*ActivatorFunc)(const G4Track&, const G4Step&);
      typedef G4ObserverCollectionT<Arg> ObserverCollection;
    };
  };  
}     
      
#endif

