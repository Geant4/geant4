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

