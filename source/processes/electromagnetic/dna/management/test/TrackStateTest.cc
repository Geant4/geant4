/*
 * TrackStateTest.cc
 *
 *  Created on: 7 avr. 2014
 *      Author: kara
 */

#include "G4TrackState.hh"
#include "globals.hh"

////////////////////////:
// Define class and State
class Class : public G4TrackStateDependent<Class>
{
public:
  Class();
  void NewTrackState();
  class State
  {
  public:
    State()
    {
      test = 0;
    }
    State(int i);
    int test;
  };
  static int ID;
  int fID;
};

int Class::ID = 0;

Class::Class()
{
  fID = Class::ID;
  Class::ID++;
}

Class::State::State(int i)
{
  test = i;
}

//LINKSTATE(Class, Class::State)

template<>
class G4TrackState<Class> : public G4TrackStateBase<Class>,
                            public Class::State
{
  typedef Class::State State;
  friend class G4TrackStateDependent<Class> ;
public:
  G4TrackState() :
      G4TrackStateBase<Class>(), Class::State()
  {
  }
  G4TrackState(int i) :
      G4TrackStateBase<Class>(), Class::State(i)
  {
  }
  virtual ~G4TrackState()
  {
  }
  virtual int GetID()
  {
    return G4TrackStateID<Class>::GetID();
  }
  static int ID()
  {
    return G4TrackStateID<Class>::GetID();
  }
protected:
};

// NewTrackState, LoadTrackState, etc... can be redefined if needed
void Class::NewTrackState()
{
  fpTrackState = StateTypeHandle(new StateType(fID));
}

// End of class and state definition
//////

void CreateState(Class& c, G4TrackStateManager& man)
{
  c.NewTrackState();
  c.SaveTrackState(man);
  c.PopTrackState();
}

int main()
{
  G4TrackStateManager man;
  Class c[10];

  for (int i = 0; i < 10; i++)
  {
    CreateState(c[i], man);
    assert(c[i].GetConcreteTrackState().get() == 0);
  }

  for (int i = 0; i < 10; i++)
  {
    G4VTrackStateHandle state1 = man.GetTrackState(&c[i]);
    G4shared_ptr<Class::State> state_converted = ConvertToConcreteTrackState<Class::State>(state1);
    assert(state_converted->test == c[i].fID);
    // G4cout << state_converted->test  << G4endl;
  }
  for (int i = 0; i < 10; i++)
  {
    c[i].LoadTrackState(man);
    Class::StateTypeHandle state = c[i].GetConcreteTrackState();
    assert(state->test == c[i].fID);
  }

  G4cout << "test OK" << G4endl;

  return 0;
}
