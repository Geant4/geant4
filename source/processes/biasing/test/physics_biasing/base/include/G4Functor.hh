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
// $Id: G4Functor.hh,v 1.1 2007-05-25 19:14:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Variation on Loki style functors:
//                         "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4FUNCTOR_HH
#define G4FUNCTOR_HH

#include "G4FunctorHandler.hh"
#include "G4FunctorImpl.hh"
#include "G4MemFunHandler.hh"
#include "G4NullFunctorHandler.hh"
#include "G4ReferenceCountedHandle.hh"

template <typename Result, typename Id, typename TList>
class G4Functor {

public:

  typedef TList ParmList;
  typedef Id Identifier;
  typedef G4FunctorImpl<Result, Identifier, TList> Impl;
  typedef typename Impl::Parm1 Parm1;
  typedef typename Impl::Parm2 Parm2;
  typedef typename Impl::Parm3 Parm3;
  typedef typename Impl::Parm4 Parm4;
  typedef typename Impl::Parm5 Parm5;
  typedef typename Impl::Parm6 Parm6;
  
  typedef G4ReferenceCountedHandle<Impl> SmartPtrImpl;

private:

  SmartPtrImpl pImpl;

public:

  typedef Result ResultType;

  // Default null handler
  G4Functor()
    :pImpl(new G4NullFunctorHandler<G4Functor>()) {}

  ~G4Functor() {}
 
  template <typename Fun> G4Functor(const Id& id, const Fun& fun) 
    :pImpl(new G4FunctorHandler<G4Functor, Fun>(id, fun)) {}
  
  G4Functor(Impl* impl) 
    :pImpl(impl) {}
  
  template <typename PtrObj, typename MemFn>
  G4Functor(const Id& id, const PtrObj&p, MemFn memFn)
    :pImpl(new G4MemFunHandler<G4Functor, PtrObj, MemFn>(id, p, memFn)) {}

  const Id& GetIdentifier() const {
    return pImpl->GetIdentifier();
  }
  /*
    //jane - to be implemented
  Id GetCompositeIdentifier() const 
  {}
  */

  Result operator()()
  {
    return (*(pImpl()))();
  }
  
  Result operator()(Parm1 p1)
  {
    return (*(pImpl()))(p1);
  }

  Result operator()(Parm1 p1, 
		    Parm2 p2)
  {
    return (*(pImpl()))(p1, p2);
  }
  
  Result operator()(Parm1 p1,
		    Parm2 p2,
		    Parm3 p3)
  {
    return (*(pImpl()))(p1, p2, p3);
  }

  Result operator()(Parm1 p1,
		    Parm2 p2,
		    Parm3 p3,
		    Parm4 p4)
  {
    return (*(pImpl()))(p1, p2, p3, p4);
  }


  Result operator()(Parm1 p1,
		    Parm2 p2,
		    Parm3 p3,
		    Parm4 p4,
		    Parm5 p5)
  {
    return (*(pImpl()))(p1, p2, p3, p4, p5);
  }

  Result operator()(Parm1 p1,
		    Parm2 p2,
		    Parm3 p3,
		    Parm4 p4,
		    Parm5 p5,
		    Parm6 p6)
  {
    return (*(pImpl()))(p1, p2, p3, p4, p5, p6);
  }
  
};

#endif
