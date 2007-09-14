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
// $Id: G4GPRTypeList.hh,v 1.5 2007-09-14 16:42:50 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. Creation. Loki style typelists:
//                       "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4GPRTYPELIST_HH
#define G4GPRTYPELIST_HH

template <class H, class T>
struct G4GPRTypeList 
{
  typedef H Head;
  typedef T Tail;
};

struct G4GPRNullType {};

struct G4GPREmptyType {};

template <typename TList, unsigned int index, typename DefaultType=G4GPRNullType> 
struct G4GPRTypeAtNonStrict
{
  typedef DefaultType Result;
};

template <typename Head, typename Tail, typename DefaultType>
struct G4GPRTypeAtNonStrict<G4GPRTypeList<Head, Tail>, 0, DefaultType>
{
  typedef Head Result;
};

template <typename Head, typename Tail, unsigned int i, typename DefaultType>
struct G4GPRTypeAtNonStrict<G4GPRTypeList<Head, Tail>, i, DefaultType>
{
  typedef typename G4GPRTypeAtNonStrict<Tail, i - 1, DefaultType>::Result  Result;
};


#define G4GPRTypeList_1(T1) G4GPRTypeList<T1, G4GPRNullType>
#define G4GPRTypeList_2(T1, T2) G4GPRTypeList<T1, G4GPRTypeList_1(T2)>
#define G4GPRTypeList_3(T1, T2, T3) G4GPRTypeList<T1, G4GPRTypeList_2(T2, T3)>
#define G4GPRTypeList_4(T1, T2, T3, T4) G4GPRTypeList<T1, G4GPRTypeList_3(T2, T3, T4)>
#define G4GPRTypeList_5(T1, T2, T3, T4, T5) G4GPRTypeList<T1, G4GPRTypeList_4(T2, T3, T4, T5)>
#define G4GPRTypeList_6(T1, T2, T3, T4, T5, T6) G4GPRTypeList<T1, G4GPRTypeList_5(T2, T3, T4, T5, T6)>
#define G4GPRTypeList_7(T1, T2, T3, T4, T5, T6, T7) G4GPRTypeList<T1, G4GPRTypeList_6(T2, T3, T4, T5, T6, T7)>
#define G4GPRTypeList_8(T1, T2, T3, T4, T5, T6, T7, T8) G4GPRTypeList<T1, G4GPRTypeList_7(T2, T3, T4, T5, T6, T7, T8)>

#endif
