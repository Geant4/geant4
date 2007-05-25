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
// $Id: G4BinderFirst.hh,v 1.1 2007-05-25 19:14:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Variation on Loki style binder 
//                         first handler :
//                         "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4BINDERFIRST_HH
#define G4BINDERFIRST_HH

template <typename A, typename B, typename C> G4Functor;

namespace {

  template <typename Functor> struct BinderFirstTraits;

  template <typename Result, typename Id, class TList>
  struct BinderFirstTraits< G4Functor<Result, Id, TList> >
  {
    typedef G4Functor<Result, Id, TList> OriginalFunctor;

    typedef typename TList::Tail ParmList;
    
    typedef typename G4TypeAtNonStrict<TList, 0>::Result OriginalParm1;
    typedef G4Functor<Result, Id, ParmList> Bound;
    typedef typename Bound::Impl Impl;
    typedef typename Bound::Identifier Identifier;
    
  };
}

template <typename Incoming>
class G4BinderFirst : public BinderFirstTraits<Incoming>::Impl {

  typedef typename Incoming::Parm2 Parm1;
  typedef typename Incoming::Parm3 Parm2;
  typedef typename Incoming::Parm4 Parm3;
  typedef typename Incoming::Parm5 Parm4;
  typedef typename Incoming::Parm6 Parm5;
  typedef  G4EmptyType Parm6;
  typedef typename BinderFirstTraits<Incoming>::Impl Base;
  typedef typename Incoming::ResultType ResultType;
  typedef typename BinderFirstTraits<Incoming>::Bound Bound;
  typedef typename BinderFirstTraits<Incoming>::Identifier Identifier;

public:

  G4BinderFirst(const Identifier& id, const Incoming& fun, Bound bound)
    :BinderFirstTraits<Incoming>::Impl(id)
    ,fFun(fun)
    ,fBound(bound) 
  {}
  
  /* jane - fixme. To be implemented
  Identifier GetCompositeIdentifier() const
  {}
  */

  ResultType operator()()
  {
    return fFun(fBound);
  }

  ResultType operator()(Parm1 p1)
  {
    return fFun(fBound, p1);
  }

  ResultType operator()(Parm1 p1, 
			Parm2 p2)
  {
    return fFun(fBound, p1, p2);
  }

  ResultType operator()(Parm1 p1, 
			Parm2 p2,
			Parm3 p3)
  {
    return fFun(fBound, p1, p2, p3);
  }

  ResultType operator()(Parm1 p1, 
			Parm2 p2,
			Parm3 p3,
			Parm4 p4)
  {
    return fFun(fBound, p1, p2, p3, p4);
  }

  ResultType operator()(Parm1 p1, 
			Parm2 p2,
			Parm3 p3,
			Parm4 p4,
			Parm5 p5)
  {
    return fFun(fBound, p1, p2, p3, p4, p5);
  }

  ResultType operator()(Parm1 p1, 
			Parm2 p2,
			Parm3 p3,
			Parm4 p4,
			Parm5 p5,
			Parm6 p6)
  {
    return fFun(fBound, p1, p2, p3, p4, p5, p6);
  }

private:
  
  Incoming fFun;
  Bound fBound;

};

template <typename Func>
G4BinderFirst<Func>*
G4BindFirst(const typename Func::Identifier& id, Func func, typename Func::Parm1 bound) 
{
  return new G4BinderFirst<Func>(id, func, bound);
}

#endif

