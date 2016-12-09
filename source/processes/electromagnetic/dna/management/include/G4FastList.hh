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
// $Id: G4FastList.hh 85244 2014-10-27 08:24:13Z gcosmo $
//
// Author: Mathieu Karamitros

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

#ifndef G4FastList_H
#define G4FastList_H

#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"
#include <G4memory.hh>
#include <vector>
#include <set>
//#include "G4ManyFastLists.hh"

template<class OBJECT>
  class G4FastList;
template<class OBJECT>
  class G4FastList_Boundary;
template<typename OBJECT>
  struct G4FastList_iterator;
template<typename OBJECT>
  struct G4FastList_const_iterator;
template<typename OBJECT>
  class G4ManyFastLists;
template<typename OBJECT>
  struct G4ManyFastLists_iterator;
template<class OBJECT>
  struct sortWatcher;


/** Comments :
 * - A track cannot belong to two different track lists
 * - Erase a given track is constant complexity
 * - This development was thought to be used together with G4IT
 */

#ifndef TYPE_WRAPPER
#define TYPE_WRAPPER
template < typename T>
struct type_wrapper
{
   typedef T type;
};
#endif

template<class LIST>
  struct _ListRef
  {
    typedef type_wrapper<LIST> traits_type;
    typedef type_wrapper<G4ManyFastLists_iterator<typename LIST::object>>
            mli_traits_type;

//#ifdef WIN32
//    friend typename traits_type::type;
//    friend typename traits_type::type::node;
//    friend typename mli_traits_type::type;
////#elif defined(__clang__)
////    friend T;
////    friend T::node;
////    friend G4ManyFastLists_iterator<traits_type::type::object>;
//#else
//    friend class traits_type::type;
//    friend class traits_type::type::node;
//    friend class mli_traits_type::type;
//#endif

    LIST* fpList;

//  protected:
    inline _ListRef(LIST* __list) :
        fpList(__list)
    {
      ;
    }
 };

/**
 * G4FastListNode is the entity actually stored
 * by the G4FastList. A G4FastListNode should
 * belong only to one list. Also, an object
 * should belong only to one list.
 */

template<class OBJECT>
  class G4FastListNode
  {
    typedef type_wrapper<OBJECT> ObjectW;
    typedef G4FastList<typename ObjectW::type> LIST;
//    typedef type_wrapper<LIST> > ListW;
    typedef type_wrapper<G4FastList<OBJECT> > ListW;
//    typedef type_wrapper<G4ManyFastLists<typename ObjectW::type> > ManyListsW;
    typedef type_wrapper<G4ManyFastLists<OBJECT> > ManyListsW;
//    typedef type_wrapper<G4ManyFastLists_iterator<typename ObjectW::type> > ManyListsIteratorW;
    typedef type_wrapper<G4ManyFastLists_iterator<OBJECT> > ManyListsIteratorW;

//#ifdef WIN32
//    friend typename ListW::type;
//    friend typename ManyListsW::type;
//    friend typename ManyListsIteratorW::type;
////#elif defined(__clang__)
////    friend T;
////    friend G4ManyFastLists_iterator<OBJECT>;
//#else
//    friend class ListW::type;
//    friend class ManyListsW::type;
//    friend struct ManyListsIteratorW::type;
//#endif

  public:
    ~G4FastListNode();

    OBJECT* GetObject()
    {
      return fpObject;
    }

    const OBJECT* GetObject() const
    {
      return fpObject;
    }

    G4FastListNode<OBJECT>* GetNext()
    {
      return fpNext;
    }
    const G4FastListNode<OBJECT>* GetNext() const
    {
      return fpNext;
    }
    G4FastListNode<OBJECT>* GetPrevious()
    {
      return fpPrevious;
    }
    const G4FastListNode<OBJECT>* GetPrevious() const
    {
      return fpPrevious;
    }
    bool IsAttached()
    {
      return fAttachedToList;
    }

  //protected:
    /** Default constructor */
    G4FastListNode(OBJECT* track = 0);

    void SetNext(G4FastListNode<OBJECT>* node)
    {
      fpNext = node;
    }
    void SetPrevious(G4FastListNode<OBJECT>* node)
    {
      fpPrevious = node;
    }
    void SetAttachedToList(bool flag)
    {
      fAttachedToList = flag;
    }

    void UnHook();

    void DetachYourSelf();
    
    bool fAttachedToList;
    G4shared_ptr<_ListRef<G4FastList<OBJECT> > > fListRef;
    OBJECT* fpObject;
    G4FastListNode<OBJECT>* fpPrevious;
    G4FastListNode<OBJECT>* fpNext;
  };

/**
 * G4FastList is used by G4TrackHolder to save
 * G4IT tracks only. Its advantage lies to a fast
 * search of a track in this list.
 */

template<class OBJECT>
  class G4FastList
  {
  protected:
    G4int fNbObjects;
//    G4FastListNode<OBJECT> * fpStart;
//    G4FastListNode<OBJECT> * fpFinish;
    G4shared_ptr<_ListRef<G4FastList<OBJECT> > > fListRef;

    G4FastListNode<OBJECT> fBoundary;
    // Must be empty and link to the last non-empty node of the list
    // and to the first non-empty node of the list (begin())
    // The iterator returned by end() is linked to this empty node

  public:
    class Watcher
    {
    public:
      enum Priority
      {
        eExtreme,
        eHigh,
        eNormal,
        eLow,
        eVeryLow
      };

      typedef G4FastList<OBJECT> list;

      Watcher()
      {
        fPriority = Priority::eVeryLow;
      }

      virtual ~Watcher()
      {
        typename std::set<G4FastList<OBJECT>*>::iterator it = fWatching.begin();
        typename std::set<G4FastList<OBJECT>*>::iterator end = fWatching.end();
        for(;it!=end;it++)
        {
          (*it)->RemoveWatcher(this);
        }
      }

      virtual G4String GetWatcherName(){
        return "";
      }

      Priority GetPriority() const{
        return fPriority;
      }

      // ===============================
      // NOTIFICATIONS
      void NotifyDeletingList(G4FastList<OBJECT>*){;}
      // used by PriorityList & ManyFastLists

      virtual void NotifyAddObject(OBJECT*, G4FastList<OBJECT>*){;}
      virtual void NotifyRemoveObject(OBJECT*, G4FastList<OBJECT>*){;}
//      void NotifyEmpty(OBJECT*, G4FastList<OBJECT>*){;}

      // ===============================

      void Watch(G4FastList<OBJECT>* fastList)
      {
        fWatching.insert(fastList);
        fastList->AddWatcher(this);
      }

      void StopWatching(G4FastList<OBJECT>* fastList, bool removeWatcher = true)
      {
        typename std::set<G4FastList<OBJECT>*>::iterator it = fWatching.find(fastList);
        if(it == fWatching.end()) return; //TODO: exception?
        fWatching.erase(it);
        if(removeWatcher) fastList->RemoveWatcher(this);
      }

    protected:
      Priority fPriority;

    private:
      std::set<G4FastList<OBJECT>*> fWatching;
    };

    template<typename WATCHER_TYPE>
    class TWatcher : public Watcher
    {
      public:
      TWatcher() : Watcher(){;}
      virtual ~TWatcher(){}
      virtual G4String GetWatcherName()
      {
        return typeid(WATCHER_TYPE).name();
      }
    };

  protected:
    typedef std::set<typename G4FastList<OBJECT>::Watcher*,
                    sortWatcher<OBJECT>> WatcherSet;
    WatcherSet   fWatchers;
    G4FastListNode<G4FastList<OBJECT> >* fpNodeInManyLists;

  public:
    typedef OBJECT object;
    typedef G4FastList_iterator<OBJECT> iterator;
    typedef G4FastList_const_iterator<OBJECT> const_iterator;
    typedef G4FastListNode<OBJECT> node;

    G4FastList();
    ~G4FastList();

    void SetListNode(G4FastListNode<G4FastList<OBJECT> >* __node)
    {
      fpNodeInManyLists = __node;
    }

    G4FastListNode<G4FastList<OBJECT> >* GetListNode()
    {
      return fpNodeInManyLists;
    }

    void AddWatcher(Watcher* watcher)
    {
      fWatchers.insert(watcher);
    }

    void RemoveWatcher(Watcher* watcher)
    {
      typename WatcherSet::iterator it = fWatchers.find(watcher);
      if(it == fWatchers.end()) return; //TODO: exception?
      fWatchers.erase(it);
    }

    inline OBJECT* back()
    {
//      if (fNbObjects != 0) return fpFinish->GetObject();
      if (fNbObjects != 0) return fBoundary.GetPrevious()->GetObject();
      else return 0;
    }

    inline G4int size() const
    {
      return fNbObjects;
    }

    inline bool empty() const;
    iterator insert(iterator /*position*/, OBJECT*);

    inline iterator begin();
    inline const_iterator begin() const;

    inline iterator end();
    inline const_iterator end() const;
    /**
     * return an iterator that contains an empty node
     * use for boundary checking only
     */

    bool Holds(const OBJECT*) const;

    inline void push_front(OBJECT* __track);
    inline void push_back(OBJECT* __track);
    OBJECT* pop_back();

    void remove(OBJECT*);

    iterator pop(OBJECT*);
    iterator pop(G4FastListNode<OBJECT>*);
    iterator pop(iterator __first, iterator __last);
    iterator erase(OBJECT*);
    /**
     * Complexity = constant
     * By storing the node inside the object, we avoid
     * searching through all the container.
     */

    iterator erase(iterator __first, iterator __last);
    /**
     * Complexity = linear in size between __first and __last
     */

    void clear();
    void transferTo(G4FastList<OBJECT>*);
    /**
     * Complexity = constant
     */

    static G4FastListNode<OBJECT>* GetNode(OBJECT*);
    static void SetNode(OBJECT* __obj, G4FastListNode<OBJECT>* __node);
    static G4FastList<OBJECT>* GetList(OBJECT*);
    static G4FastList<OBJECT>* GetList(G4FastListNode<OBJECT>* __trackListNode);
    static void Pop(OBJECT*);

  protected:
    G4FastListNode<OBJECT>* CreateNode(OBJECT*);
    static G4FastListNode<OBJECT>* __GetNode(OBJECT*);
    G4FastListNode<OBJECT>* Flag(OBJECT*);
    G4FastListNode<OBJECT>* Unflag(OBJECT*);
    void Unflag(G4FastListNode<OBJECT>* __trackListNode);
    void CheckFlag(G4FastListNode<OBJECT>*);
    void DeleteObject(OBJECT*);

    void Hook(G4FastListNode<OBJECT>* /*position*/,
              G4FastListNode<OBJECT>* /*toHook*/);
    void Unhook(G4FastListNode<OBJECT>*);
    G4FastListNode<OBJECT>* EraseListNode(OBJECT*);

  private:
    G4FastList(const G4FastList<OBJECT>& other);
    G4FastList<OBJECT> & operator=(const G4FastList<OBJECT> &right);
    G4int operator==(const G4FastList<OBJECT> &right) const;
    G4int operator!=(const G4FastList<OBJECT> &right) const;
  };


template<class OBJECT>
  struct sortWatcher
  {
    bool operator()(const typename G4FastList<OBJECT>::Watcher* left,
                    const typename G4FastList<OBJECT>::Watcher* right) const
    {
      if(left && right)
      {
        if(left->GetPriority() != right->GetPriority())
        {
          return left->GetPriority() < right->GetPriority();
        }
        return left < right;
      }
      return false;
    }
  };


/**
 * G4FastList_iterator enables to go through
 * the tracks contained by a list.
 */

template<typename OBJECT>
  struct G4FastList_iterator
  {
//    friend class G4FastList<OBJECT>;
    typedef G4FastList_iterator<OBJECT> _Self;
    typedef G4FastListNode<OBJECT> _Node;

    G4FastList_iterator() :
        fpNode(0)
    {
    }

    explicit G4FastList_iterator(_Node* __x) :
        fpNode(__x)
    {
    }

    G4FastList_iterator(const G4FastList_iterator& right) :
        fpNode(right.fpNode)
    {
    }

    _Node* GetNode()
    {
      return fpNode;
    }

    const _Node* GetNode() const
    {
      return fpNode;
    }

    OBJECT*
    operator*();

    const OBJECT*
    operator*() const;

    OBJECT*
    operator->();

    const OBJECT*
    operator->() const;

    _Self&
    operator++()
    {
      fpNode = fpNode->GetNext();
      return *this;
    }

    _Self operator++(int)
    {
      _Self __tmp = *this;
      fpNode = fpNode->GetNext();
      return __tmp;
    }

    _Self&
    operator--()
    {
      fpNode = fpNode->GetPrevious();
      return *this;
    }

    _Self operator--(int)
    {
      _Self __tmp = *this;
      fpNode = fpNode->GetPrevious();
      return __tmp;
    }

    bool operator==(const _Self& __x) const
    {
      return (fpNode == __x.fpNode);
    }

    bool operator!=(const _Self& __x) const
    {
      return (fpNode != __x.fpNode);
    }

//  private:
    // The only member points to the G4FastList_iterator element.
    _Node* fpNode;
  };

/**
 * G4FastList_iterator enables to go through
 * the tracks contained by a list.
 */

template<typename OBJECT>
  struct G4FastList_const_iterator
  {
//    friend class G4FastList<OBJECT>;
    typedef G4FastList_const_iterator<OBJECT> _Self;
    typedef G4FastListNode<OBJECT> _Node;

    G4FastList_const_iterator() :
        fpNode(0)
    {
    }

    explicit G4FastList_const_iterator(const _Node* __x) :
        fpNode(__x)
    {
    }

    G4FastList_const_iterator(const G4FastList_const_iterator& right) :
        fpNode(right.fpNode)
    {
    }

    G4FastList_const_iterator(const G4FastList_iterator<OBJECT>& right) :
        fpNode(right.GetNode())
    {
    }

    const OBJECT*
    operator*() const
    {
      if(fpNode == 0) return 0;
      return fpNode->GetObject();
    }

    const OBJECT*
    operator->() const
    {
      if(fpNode == 0) return 0;
      return fpNode->GetObject();
    }

    _Self&
    operator++()
    {
      fpNode = fpNode->GetNext();
      return *this;
    }

    _Self operator++(int)
    {
      _Self __tmp = *this;
      fpNode = fpNode->GetNext();
      return __tmp;
    }

    _Self&
    operator--()
    {
      fpNode = fpNode->GetPrevious();
      return *this;
    }

    _Self operator--(int)
    {
      _Self __tmp = *this;
      fpNode = fpNode->GetPrevious();
      return __tmp;
    }

    bool operator==(const _Self& __x) const
    {
      return (fpNode == __x.fpNode);
    }

    bool operator!=(const _Self& __x) const
    {
      return (fpNode != __x.fpNode);
    }

//  private:
    // The only member points to the G4FastList_iterator element.
    const _Node* fpNode;
  };

#include "G4FastList.icc"

#endif // G4FastList_H
