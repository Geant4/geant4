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
// $Id: G4TypeList.hh,v 1.1 2007-05-25 19:14:37 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Loki style typelists:
//                         "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4TYPELIST_HH
#define G4TYPELIST_HH

template <class H, class T>
struct G4TypeList 
{
  typedef H Head;
  typedef T Tail;
};

struct G4NullType {};

struct G4EmptyType {};

template <typename TList, unsigned int index, typename DefaultType=G4NullType> 
struct G4TypeAtNonStrict
{
  typedef DefaultType Result;
};

template <typename Head, typename Tail, typename DefaultType>
struct G4TypeAtNonStrict<G4TypeList<Head, Tail>, 0, DefaultType>
{
  typedef Head Result;
};

template <typename Head, typename Tail, unsigned int i, typename DefaultType>
struct G4TypeAtNonStrict<G4TypeList<Head, Tail>, i, DefaultType>
{
  typedef typename G4TypeAtNonStrict<Tail, i - 1, DefaultType>::Result  Result;
};

#define G4TypeList_1(T1) G4TypeList<T1, G4NullType>
#define G4TypeList_2(T1, T2) G4TypeList<T1, G4TypeList_1(T2)>
#define G4TypeList_3(T1, T2, T3) G4TypeList<T1, G4TypeList_2(T2, T3)>
#define G4TypeList_4(T1, T2, T3, T4) G4TypeList<T1, G4TypeList_3(T2, T3, T4)>
#define G4TypeList_5(T1, T2, T3, T4, T5) G4TypeList<T1, G4TypeList_4(T2, T3, T4, T5)>
#define G4TypeList_6(T1, T2, T3, T4, T5, T6) G4TypeList<T1, G4TypeList_5(T2, T3, T4, T5, T6)>

#endif
