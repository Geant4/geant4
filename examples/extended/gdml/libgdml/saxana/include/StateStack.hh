//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: StateStack.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef STATE_STACK_H
#define STATE_STACK_H 1

//#include "RCObjectHandle.hh"
#include "ProcessingState.hh"

#include <stack>

class StateStack
{
public:
  //typedef RCObjectHandle<ProcessingState> State;
  //typedef std::stack<State>               Stack;
  typedef ProcessingState*      State;
  typedef std::stack<State>     Stack;

public:
  bool Empty() const;
  unsigned int Size() const;
  StateStack::State& Top();
  const StateStack::State& Top() const;
  void Push( const StateStack::State& s );
  //void Push( ProcessingState* s );
  void Pop();
  
private:
  Stack  fStack;
};

inline bool StateStack::Empty() const { return fStack.empty(); }
inline unsigned int StateStack::Size() const { return fStack.size(); }
inline StateStack::State& StateStack::Top() { return fStack.top(); }
inline const StateStack::State& StateStack::Top() const { return fStack.top(); }
inline void StateStack::Push( const StateStack::State& s ) { fStack.push( s ); }
//inline void StateStack::Push( ProcessingState* s )
//{
//  State st = s;
//  fStack.push( st );
//}
inline void StateStack::Pop() { fStack.pop(); }


#endif // STATE_STACK_H

