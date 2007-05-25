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
// $Id: G4FunctorImpl.hh,v 1.1 2007-05-25 19:14:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Variation on Loki style functor handler:
//                         "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4FUNCTORIMPL_HH
#define G4FUNCTORIMPL_HH

#include "G4TypeList.hh"
#include "G4TypeTraits.hh"

namespace {

  template <typename Id>
  class FunctorImplBase {
  public:
    FunctorImplBase(const Id& id): fId(id) {}
    
    typedef Id Identifier;
    typedef G4EmptyType Parm1;
    typedef G4EmptyType Parm2;
    typedef G4EmptyType Parm3;
    typedef G4EmptyType Parm4;
    typedef G4EmptyType Parm5;
    typedef G4EmptyType Parm6;

    virtual const Id& GetIdentifier() const {return fId;}

    /* jane fixme - to be implemented
    virtual Id GetCompositeIdentifier() const 
    {}
    */
    
  private:
    Id fId;
  };
}

template <typename Result, typename Identifier, typename TList>
struct G4FunctorImpl;

template <typename Result, typename Identifier>
struct G4FunctorImpl<Result, Identifier, G4NullType> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}

  virtual ~G4FunctorImpl() {}
  virtual Result operator()() = 0;
  
};

template <typename Result, typename Identifier, typename P1>
struct G4FunctorImpl<Result, Identifier, G4TypeList_1(P1)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  virtual Result operator()(Parm1) = 0;

};

template <typename Result, typename Identifier, typename P1, typename P2>
struct G4FunctorImpl<Result, Identifier, G4TypeList_2(P1, P2)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  typedef typename G4TypeTraits<P2>::ParameterType Parm2;

  virtual Result operator()(Parm1, 
			    Parm2) = 0;
  
};

template <typename Result, typename Identifier, typename P1, typename P2, typename P3>
struct G4FunctorImpl<Result, Identifier, G4TypeList_3(P1, P2, P3)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  typedef typename G4TypeTraits<P2>::ParameterType Parm2;
  typedef typename G4TypeTraits<P3>::ParameterType Parm3;

  virtual Result operator()(Parm1, 
			    Parm2, 
			    Parm3) = 0;
  
};

template <typename Result, typename Identifier, typename P1, typename P2, typename P3, typename P4>
struct G4FunctorImpl<Result, Identifier, G4TypeList_4(P1, P2, P3, P4)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  typedef typename G4TypeTraits<P2>::ParameterType Parm2;
  typedef typename G4TypeTraits<P3>::ParameterType Parm3;
  typedef typename G4TypeTraits<P4>::ParameterType Parm4;

  virtual Result operator()(Parm1, 
			    Parm2, 
			    Parm3,
			    Parm4) = 0;
  
};

template <typename Result, typename Identifier, typename P1, typename P2, typename P3, typename P4, typename P5>
struct G4FunctorImpl<Result, Identifier, G4TypeList_5(P1, P2, P3, P4, P5)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  typedef typename G4TypeTraits<P2>::ParameterType Parm2;
  typedef typename G4TypeTraits<P3>::ParameterType Parm3;
  typedef typename G4TypeTraits<P4>::ParameterType Parm4;
  typedef typename G4TypeTraits<P5>::ParameterType Parm5;

  virtual Result operator()(Parm1, 
			    Parm2, 
			    Parm3,
			    Parm4,
			    Parm5) = 0;
  
};

template <typename Result, typename Identifier, typename P1, typename P2, typename P3, typename P4, typename P5, typename P6>
struct G4FunctorImpl<Result, Identifier, G4TypeList_6(P1, P2, P3, P4, P5, P6)> : public FunctorImplBase<Identifier> {

public:
  
  G4FunctorImpl(Identifier id): FunctorImplBase<Identifier>(id){}
  virtual ~G4FunctorImpl() {}

  typedef typename G4TypeTraits<P1>::ParameterType Parm1;
  typedef typename G4TypeTraits<P2>::ParameterType Parm2;
  typedef typename G4TypeTraits<P3>::ParameterType Parm3;
  typedef typename G4TypeTraits<P4>::ParameterType Parm4;
  typedef typename G4TypeTraits<P5>::ParameterType Parm5;
  typedef typename G4TypeTraits<P6>::ParameterType Parm6;

  virtual Result operator()(Parm1, 
			    Parm2, 
			    Parm3,
			    Parm4,
			    Parm5,
			    Parm6) = 0;
  
};

#endif
