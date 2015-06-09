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
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4TRACKSTATE_HH_
#define G4TRACKSTATE_HH_

#include "G4memory.hh"
#include <map>

class G4VTrackStateID
{
public:
  virtual int GetStateID() = 0;
protected:
  G4VTrackStateID()
  {
    ;
  }
  virtual ~G4VTrackStateID()
  {
    ;
  }

  static int fgLastID;
};

template<class T>
  class G4TrackStateID : public G4VTrackStateID
  {
  public:
    static int GetID()
    {
      return fID;
    }

  private:
    static int Create()
    {
      fID = fgLastID;
      fgLastID++;
      return fID;
    }

    G4TrackStateID()
    {
    }
    ~G4TrackStateID()
    {
      ;
    }
    static int fID;
  };

template<class T>
  int G4TrackStateID<T>::fID = G4TrackStateID<T>::Create();

class G4VTrackState
{
public:
  G4VTrackState()
  {
  }
  virtual ~G4VTrackState()
  {
  }
  virtual int GetID() = 0;
};

typedef G4shared_ptr<G4VTrackState> G4VTrackStateHandle;

//template<class TrackStateDependent>
//class RegisterState : public TrackStateDependent
//{
//public:
//  typedef typename TrackStateDependent::State TrackState;
//};

#ifndef TYPE_WRAPPER
#define TYPE_WRAPPER
template < typename T>
struct type_wrapper
{
   typedef T type;
};
#endif

//!
template<class T>
  class G4TrackStateDependent;

template<class T>
  class G4TrackStateBase : public G4VTrackState
  {
    //friend class type_wrapper<T>::type; // works with gcc
    typedef type_wrapper<T> traits_type;
    /*
     #ifdef WIN32
     friend typename traits_type::type;
     #elif defined(__clang__)
     friend T;
     #else
     //        friend class traits_type::type;
     friend class type_wrapper<T>::type;
     #endif
     */
    friend class G4TrackStateDependent<T> ; //!

  public:
    virtual ~G4TrackStateBase()
    {
    }

    virtual int GetID()
    {
      return G4TrackStateID<T>::GetID();
    }

    static int ID()
    {
      return G4TrackStateID<T>::GetID();
    }

    G4TrackStateBase() :
        G4VTrackState()
    {
    }

  protected:

  };

template<class T>
  class G4TrackState : public G4TrackStateBase<T>
  {
    //  friend class type_wrapper<T>::type; // works with gcc

    typedef type_wrapper<T> traits_type;
    /*
     #ifdef WIN32
     friend typename traits_type::type;
     #elif defined(__clang__)
     friend T;
     #else
     friend class type_wrapper<T>::type;
     //friend class traits_type::type;
     #endif
     */
    friend class G4TrackStateDependent<T> ; //!
  public:
    virtual ~G4TrackState()
    {
    }

    virtual int GetID()
    {
      return G4TrackStateID<T>::GetID();
    }

    static int ID()
    {
      return G4TrackStateID<T>::GetID();
    }

    G4TrackState() :
        G4TrackStateBase<T>()
    {
    }

  protected:
  };

class G4TrackStateManager
{
  std::map<int, G4VTrackStateHandle> fTrackStates;
  std::map<void*, G4VTrackStateHandle> fMultipleTrackStates;

public:

  void SetTrackState(void* adress, G4VTrackStateHandle state)
  {
    fMultipleTrackStates[adress] = state;
  }

  G4VTrackStateHandle GetTrackState(void* adress) const
  {
    std::map<void*, G4VTrackStateHandle>::const_iterator it =
        fMultipleTrackStates.find(adress);
    if (it == fMultipleTrackStates.end())
    {
      return G4VTrackStateHandle();
    }
    return it->second;
  }

  template<class T>
    G4VTrackStateHandle GetTrackState(T* adress) const
    {
      std::map<void*, G4VTrackStateHandle>::const_iterator it =
          fMultipleTrackStates.find((void*) adress);
      if (it == fMultipleTrackStates.end())
      {
        return G4VTrackStateHandle();
      }
      return it->second;
    }
  void SetTrackState(G4VTrackStateHandle state)
  {
    fTrackStates[state->GetID()] = state;
  }
  template<typename T>
    G4VTrackStateHandle GetTrackState() const
    {
      std::map<int, G4VTrackStateHandle>::const_iterator it = fTrackStates.find(
          G4TrackStateID<T>::GetID());
      if (it == fTrackStates.end())
      {
        return G4VTrackStateHandle();
      }
      return it->second;
    }
};

class G4VTrackStateDependent
{
public:
  G4VTrackStateDependent()
  {
    ;
  }
  virtual ~G4VTrackStateDependent()
  {
    ;
  }

  virtual void NewTrackState() = 0;
  virtual void LoadTrackState(G4TrackStateManager&) = 0;
  virtual void SaveTrackState(G4TrackStateManager&) = 0;
  virtual G4VTrackStateHandle GetTrackState() const = 0;
  virtual G4VTrackStateHandle PopTrackState() = 0;
  virtual void ResetTrackState() = 0;
};

#define G4TrackStateHandle(T) G4shared_ptr<G4TrackState<T> >

template<class OriginalType>
  G4shared_ptr<G4VTrackState> ConvertToAbstractTrackState(G4shared_ptr<G4TrackState<OriginalType> > state)
  {

    G4shared_ptr<G4VTrackState> output = G4dynamic_pointer_cast<G4VTrackState>(state);
    return output;
  }

template<class FinalType>
  G4shared_ptr<G4TrackState<FinalType> > ConvertToConcreteTrackState(G4VTrackStateHandle state)
  {

    G4shared_ptr<G4TrackState<FinalType> > output = G4dynamic_pointer_cast<G4TrackState<FinalType> >(state);
    return output;
  }

//!
template<class T>
  class G4TrackStateDependent : public G4VTrackStateDependent
  {
  public:
    typedef T ClassType;
    typedef G4TrackState<T> StateType;
    typedef G4shared_ptr<StateType> StateTypeHandle;

    virtual ~G4TrackStateDependent()
    { ;}
    virtual void SetTrackState(G4shared_ptr<StateType> state)
    {
      fpTrackState = state;
    }

    virtual G4VTrackStateHandle PopTrackState()
    {
      G4VTrackStateHandle output = G4dynamic_pointer_cast<G4VTrackState>(fpTrackState);
      fpTrackState.reset();
      return output;
    }

    virtual G4VTrackStateHandle GetTrackState() const
    {
      G4VTrackStateHandle output = G4dynamic_pointer_cast<G4VTrackState>(fpTrackState);
      return output;
    }

    virtual StateTypeHandle GetConcreteTrackState() const
    {
      return fpTrackState;
    }

    virtual void LoadTrackState(G4TrackStateManager& manager)
    {
      fpTrackState = ConvertToConcreteTrackState<ClassType>(manager.GetTrackState(this));
      if(fpTrackState == 0)
      {
        NewTrackState();
        SaveTrackState(manager);
      }
    }

    virtual void SaveTrackState(G4TrackStateManager& manager)
    {
      manager.SetTrackState(this, ConvertToAbstractTrackState(fpTrackState));
    }

    virtual void NewTrackState()
    {
      fpTrackState = StateTypeHandle(new StateType());
    }

    virtual StateTypeHandle CreateTrackState() const
    {
      return StateTypeHandle(new StateType());
    }

    virtual void ResetTrackState()
    {
      fpTrackState.reset();
    }

  protected:
    G4TrackStateDependent()
    : G4VTrackStateDependent()
    { ;}

    StateTypeHandle fpTrackState;
  };

#define RegisterTrackState(CLASS,STATE) \
    template<> \
    class G4TrackState<CLASS> : public G4TrackStateBase<CLASS>, public STATE \
    { \
  typedef STATE State; \
  friend class G4TrackStateDependent<CLASS>;  \
    public: \
    G4TrackState() : G4TrackStateBase<CLASS>(),STATE(){}\
    virtual ~G4TrackState(){}\
    virtual int GetID()\
    {\
      return  G4TrackStateID<CLASS>::GetID();\
    }\
    static int ID()\
    {\
      return  G4TrackStateID<CLASS>::GetID();\
    }\
    protected:\
    };\


#endif /* G4TRACKSTATE_HH_ */
