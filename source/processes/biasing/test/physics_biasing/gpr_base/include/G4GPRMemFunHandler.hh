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
// $Id: G4GPRMemFunHandler.hh,v 1.2 2007-09-06 22:07:04 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. Creation. Variation on Loki style 
//                       member function functor handler:
//                       "Modern C++ Design, Andrei Alexandrescu"
#ifndef G4GPRMEMFUNHANDLER_HH
#define G4GPRMEMFUNHANDLER_HH

template <typename ParentFunctor, typename PointerToObj, 
	  typename PointerToMemFn>
class G4GPRMemFunHandler : public ParentFunctor::Impl {
  
  PointerToObj fPObj;
  PointerToMemFn fPMemFn;
  
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
  
  G4GPRMemFunHandler(const Identifier& id, const PointerToObj& pObj, PointerToMemFn pMemFn)
    :ParentFunctor::Impl(id)
    ,fPObj(pObj)
    ,fPMemFn(pMemFn)
  {}

  ResultType operator()() 
  {
    return ((*fPObj).*fPMemFn)();
  }

  ResultType operator()(Parm1 p1) 
  {
    return ((*fPObj).*fPMemFn)(p1);
  }

  ResultType operator()(Parm1 p1, Parm2 p2) 
  {
    return ((*fPObj).*fPMemFn)(p1, p2);
  }

  ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3) 
  {
    return ((*fPObj).*fPMemFn)(p1, p2, p3);
  }

  ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4) 
  {
    return ((*fPObj).*fPMemFn)(p1, p2, p3, p4);
  }

  ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5) 
  {
    return ((*fPObj).*fPMemFn)(p1, p2, p3, p4, p5);
  }

  ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5, Parm6 p6) 
  {
    return ((*fPObj).*fPMemFn)(p1, p2, p3, p4, p5, p6);
  }

};
#endif
