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
// $Id: G4GPRTriggerTypes.hh,v 1.1 2007-09-10 22:05:01 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// J. Tinslay, May 2007. Creation - scope definitions.
//
#ifndef G4GPRTRIGGERTYPES_HH
#define G4GPRTRIGGERTYPES_HH
   
#include "G4String.hh"
#include "G4GPRObserverCollectionT.hh"
#include "G4GPRTypeList.hh"
class G4Track;
class G4Step;
    
namespace G4GPRTriggerTypes {
    
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
    typedef G4GPRFunctor<G4bool, G4String, Args> ActivatorDef;
    typedef G4GPRObserverCollectionT<Args> ObserverCollection;
  };
  struct Global {
    typedef G4GPRNullType Arg;
    typedef G4GPRFunctor<G4bool, G4String, Arg> ActivatorDef;  
    typedef G4bool (*ActivatorFunc)();
    typedef G4GPRObserverCollectionT<Arg> ObserverCollection;
  };

  struct Tracking {
    struct StartTracking {
      typedef G4GPRTypeList_1(G4Track*) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;
      typedef G4bool (*TriggerFunc)(G4Track*);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(G4Track*);
      };
    };

    struct EndTracking {
      typedef G4GPRTypeList_1(G4Track) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> ActivatorDef;
      typedef G4bool (*TriggerFunc)(const G4Track&);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&);
      };
    };
  };
     
  struct Stepping {
    struct StartStep {
      typedef G4GPRTypeList_2(G4Track, G4Step) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;

      typedef G4bool (*TriggerFunc)(const G4Track&, const G4Step&);
      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&, const G4Step&);
      };
    };
      
    struct EndStep {
      typedef G4GPRTypeList_2(G4Track, G4Step) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;

      typedef G4bool (*TriggerFunc)(const G4Track&, const G4Step&);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&, const G4Step&);
      };
    };
  };
  struct Geometry {
    struct StartBoundary {
      typedef G4GPRTypeList_2(G4Track, G4Step) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;
      typedef G4bool (*TriggerFunc)(const G4Track&, const G4Step&);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&, const G4Step&);
      };
    };
      
    struct EndBoundary {
      typedef G4GPRTypeList_2(G4Track, G4Step) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;
      typedef G4bool (*TriggerFunc)(const G4Track&, const G4Step&);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&, const G4Step&);
      };
    };

    struct NewRegion {
      typedef G4GPRTypeList_1(G4Track) Arg;
      typedef G4GPRFunctor<G4bool, G4String, Arg> TriggerWrapper;
      typedef G4bool (*TriggerFunc)(const G4Track&);

      template <typename T>
      struct TriggerMfn {
	typedef G4bool (T::*PtrToMfn)(const G4Track&);
      };
    };

  };  
}     
      
#endif

