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
// $Id: G4FunctorHandler.hh,v 1.1 2007-05-25 19:14:36 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Variation on Loki style functor handler:
//                         "Modern C++ Design, Andrei Alexandrescu"
#ifndef G4FUNCTORHANDLER_HH
#define G4FUNCTORHANDLER_HH

#include "G4FunctorImpl.hh"

template <typename ParentFunctor, typename Fun>
class G4FunctorHandler : public ParentFunctor::Impl {

public:

  typedef typename ParentFunctor::Impl Base;
  typedef typename ParentFunctor::ResultType ResultType;

  typedef typename Base::Identifier Identifier;
  typedef typename Base::Parm1 Parm1;
  typedef typename Base::Parm2 Parm2;
  typedef typename Base::Parm3 Parm3;
  typedef typename Base::Parm4 Parm4;
  typedef typename Base::Parm5 Parm5;
  typedef typename Base::Parm6 Parm6;

  G4FunctorHandler(const Identifier& id, const Fun& fun) 
    :ParentFunctor::Impl(id)
    ,fFun(fun) {}

  /* jane fixme - to be implemented
  Identifier GetCompositeIdentifier() const
  {}
  */

  ResultType operator()() 
  {
    return fFun();
  }

  ResultType operator()(Parm1 p1)
  {
    return fFun(p1);
  }

  ResultType operator()(Parm1 p1,
			Parm2 p2)
  {
    return fFun(p1, p2);
  }

  ResultType operator()(Parm1 p1,
			Parm2 p2,
			Parm3 p3)
  {
    return fFun(p1, p2, p3);
  }

  ResultType operator()(Parm1 p1,
			Parm2 p2,
			Parm3 p3,
			Parm4 p4)
  {
    return fFun(p1, p2, p3, p4);
  }

  ResultType operator()(Parm1 p1,
			Parm2 p2,
			Parm3 p3,
			Parm4 p4,
			Parm5 p5)
  {
    return fFun(p1, p2, p3, p4, p5);
  }

  ResultType operator()(Parm1 p1,
			Parm2 p2,
			Parm3 p3,
			Parm4 p4,
			Parm5 p5,
			Parm6 p6)
  {
    return fFun(p1, p2, p3, p4, p5, p6);
  }



private:

  Fun fFun;

};

#endif
