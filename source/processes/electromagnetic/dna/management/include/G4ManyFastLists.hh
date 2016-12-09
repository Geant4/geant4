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
/*
 * G4ManyFastLists.hh
 *
 *  Created on: 17 nov. 2014
 *      Author: kara
 */

#ifndef G4MANYFASTLISTS_HH_
#define G4MANYFASTLISTS_HH_

#include "G4FastList.hh"
#include <set>

template<class OBJECT>
  struct G4ManyFastLists_iterator;

/*
 * Roll over many list as if it was one
 */
template<class OBJECT>
  class G4ManyFastLists : public G4FastList<OBJECT>::Watcher
  {
  protected:
    typedef G4FastList<G4FastList<OBJECT> > ManyLists;
    ManyLists fAssociatedLists;
    // TODO use "marked list" insted of vector

    typedef std::set<typename G4FastList<OBJECT>::Watcher*,
        sortWatcher<OBJECT>> WatcherSet;
    WatcherSet* fMainListWatchers;

  public:
    typedef G4ManyFastLists_iterator<OBJECT> iterator;

    G4ManyFastLists() : G4FastList<OBJECT>::Watcher(),
        fAssociatedLists(), fMainListWatchers(0)
    {
    }

    virtual ~G4ManyFastLists(){;}

    virtual void NotifyDeletingList(G4FastList<OBJECT>* __list)
    {
      fAssociatedLists.pop(__list);
    }

    void AddGlobalWatcher(typename G4FastList<OBJECT>::Watcher* watcher)
    {
      if(fMainListWatchers == 0)
      {
        fMainListWatchers = new WatcherSet();
      }

      G4cout << watcher->GetWatcherName() << G4endl;

      fMainListWatchers->insert(watcher);

      typename ManyLists::iterator it = fAssociatedLists.begin();
      typename ManyLists::iterator _end = fAssociatedLists.end();

      for(;it != _end ;++it)
      {
        watcher->Watch(*it);
//        (*it)->AddWatcher(watcher);
//        (*it)->AddWatcher(watcher);
      }
    }

    inline void Add(G4FastList<OBJECT>* __list)
    {
      if (__list == 0) return;
      fAssociatedLists.push_back(__list); // TODO use the table doubling tech
      //__list->AddWatcher(this);
      this->Watch(__list);

      if(fMainListWatchers == 0) return;

      typename WatcherSet::iterator it_watcher = fMainListWatchers->begin();
      typename WatcherSet::iterator end_watcher = fMainListWatchers->end();

//      G4cout << "G4ManyFastLists::Add -- N watchers ="
//             << fMainListWatchers->size()
//             << G4endl;

      for(;it_watcher != end_watcher ;++it_watcher)
      {
//        G4cout << " *** *** *** WATCH --- "
//                           << (*it_watcher)->GetWatcherName()
//                           << G4endl;
        (*it_watcher)->Watch(__list);
      }

      if(__list->empty() == false)
      {
        it_watcher = fMainListWatchers->begin();

        for(;it_watcher != end_watcher ;++it_watcher)
        {
          typename G4FastList<OBJECT>::iterator it_obj = __list->begin();
          for(;it_obj != __list->end() ;++it_obj)
          {
//            G4cout << " *** *** *** NOTIFY ADD OBJ --- "
//                   << (*it_watcher)->GetWatcherName()
//                   << G4endl;

            (*it_watcher)->NotifyAddObject(*it_obj,__list);
          }
        }
      }
//      else
//      {
//        G4cout << "__list->empty() == true" << G4endl;
//      }

      /*
      typename ManyLists::const_iterator __it = fAssociatedLists
          .begin();
      typename ManyLists::const_iterator __end = fAssociatedLists
          .end();
      for (; __it != __end; __it++)
      {
        assert(*__it);
      }
      */
    }

    inline void Remove(G4FastList<OBJECT>* __list)
    {
      if (__list == 0) return;
      fAssociatedLists.pop(__list); // TODO use the table doubling tech
      __list->RemoveWatcher(this);
      this->StopWatching(__list);

      typename WatcherSet::iterator it = fMainListWatchers->begin();
      typename WatcherSet::iterator _end = fMainListWatchers->end();

      for(;it != _end ;++it)
      {
        (*it)->StopWatching(__list);
      }

//      typename ManyLists::node* __node = __list->GetListNode();
//      if(__node)
//      {
//        __list->SetListNode(0);
//        delete __node;
//      }
    }

    inline bool Holds(OBJECT* __track) const
    {
      typename ManyLists::const_iterator __it = fAssociatedLists.begin();
      typename ManyLists::const_iterator __end = fAssociatedLists.end();
      for (; __it != __end; __it++)
        if ((*__it)->Holds(__track)) return true;
      return false;
    }

    inline size_t size() const
    {
      size_t __size(0);
      typename ManyLists::const_iterator __it = fAssociatedLists
          .begin();
      typename ManyLists::const_iterator __end = fAssociatedLists
          .end();
      for (; __it != __end; __it++)
      {
        __size += (*__it)->size();
      }
      return __size;
    }

    inline void RemoveLists()
    {
      typename ManyLists::iterator __it = fAssociatedLists.begin();
      typename ManyLists::iterator __end = fAssociatedLists.end();
      for (; __it != __end; __it++)
      {
        if (*__it)
        {
          (*__it)->clear();
          typename ManyLists::iterator next = __it;
          next++;
          Remove(*__it);
          typename ManyLists::node* __node = __it.GetNode();
          if(__node)
          {
            __node->GetObject()->SetListNode(0);
            delete __node;
          }
//          delete (*__it);

          __it = next;
        }
      }
      fAssociatedLists.clear();
    }

    inline void ClearLists()
    {
      typename ManyLists::iterator __it = fAssociatedLists.begin();
      typename ManyLists::iterator __end = fAssociatedLists.end();
      for (; __it != __end; __it++)
        if (*__it) (*__it)->clear();
    }

    inline iterator begin();
    inline iterator end();

    void pop(OBJECT*);
  };

template<class OBJECT>
  struct G4ManyFastLists_iterator
  {
//    friend class G4ManyFastLists<OBJECT>;
    typedef G4FastList<G4FastList<OBJECT> > ManyLists;

    typedef G4ManyFastLists_iterator _Self;
    typedef G4FastListNode<OBJECT> _Node;

    G4FastList_iterator<OBJECT> fIterator;
    typename ManyLists::iterator fCurrentListIt;
    ManyLists* fLists;

  private:
    G4ManyFastLists_iterator() :
        fIterator(), fLists(0)
    {
    }

  public:

    explicit G4ManyFastLists_iterator(G4FastList_iterator<OBJECT> __x,
                                      typename ManyLists::iterator __it,
                                      ManyLists* __lists) :
        fIterator(__x), fCurrentListIt(__it), fLists(__lists)
    {
    }

    G4ManyFastLists_iterator(const G4ManyFastLists_iterator& __x) :
        fIterator(__x.fIterator),
        fCurrentListIt(__x.fCurrentListIt),
        fLists(__x.fLists)
    {
    }

    _Node* GetNode()
    {
      return fIterator.GetNode();
    }

    G4FastList<OBJECT>* GetTrackList()
    {
      return *fCurrentListIt;
    }

    OBJECT* operator*()
    {
      return *fIterator;
    }
    const OBJECT* operator*() const
    {
      return *fIterator;
    }
    OBJECT* operator->()
    {
      return *fIterator;
    }
    const OBJECT* operator->() const
    {
      return *fIterator;
    }

    _Self UpdateToNextValidList();
    _Self& operator++();

    _Self operator++(int)
    {
      return operator++();
    }

    _Self&
    operator--()
    {
      if (fLists->empty())
      {
        fIterator = G4FastList_iterator<OBJECT>();
        return *this;
      }
      if (fCurrentListIt == fLists->begin())
      {
        if (fIterator == (*fCurrentListIt)->begin())
        {
          fIterator = G4FastList_iterator<OBJECT>();
          return *this;
        }
      }

      if (fCurrentListIt == fLists->end())
      {
        fCurrentListIt--;
        fIterator = (*fCurrentListIt)->end();
      }
      else if (fIterator == (*fCurrentListIt)->begin())
      {
        fCurrentListIt--;
        fIterator = (*fCurrentListIt)->end();
      }

      fIterator--;

      while (((*fCurrentListIt)->empty() || fIterator.GetNode() == 0
              || fIterator.GetNode()->GetObject() == 0)
             && fCurrentListIt != fLists->begin())
      {
        fIterator = (*fCurrentListIt)->begin();
        fCurrentListIt--;
        fIterator = (*fCurrentListIt)->end();
        fIterator--;
      }

      if (fIterator.GetNode() == 0 && fCurrentListIt == fLists->begin())
      {
        fIterator = G4FastList_iterator<OBJECT>();
        return *this;
      }

      return *this;
    }

    _Self operator--(int)
    {
      return operator--();
    }

    bool operator==(const _Self& __x) const
    {
      return (fIterator == __x.fIterator && fCurrentListIt == __x.fCurrentListIt);
    } // Fast check

    bool operator!=(const _Self& __x) const
    {
      return !(this->operator ==(__x));
    }

  protected:
    void HasReachedEnd()
    {
      if (fLists->empty() == false)
      {
        fIterator = (*(fLists->end()--))->end();
      }
      else
      {
        fIterator = G4FastList_iterator<OBJECT>();
      }
    }
  };

template<class OBJECT>
  typename G4ManyFastLists<OBJECT>::iterator G4ManyFastLists<OBJECT>::begin()
  {
    if (fAssociatedLists.empty())
    {
      return G4ManyFastLists_iterator<OBJECT>(G4FastList_iterator<OBJECT>(),
                                             fAssociatedLists.end(),
                                             &fAssociatedLists);
    }

    typename G4FastList<OBJECT>::iterator trackList_it;
    int i = 0;

    typename ManyLists::iterator it = fAssociatedLists.begin();
    typename ManyLists::iterator _end = fAssociatedLists.end();

    while (it != _end)
    {
      if (*it && (*it)->empty() == false)
      {
        trackList_it = (*it)->begin();
        break;
      }
      i++;
      it++;
    };

    if (i == fAssociatedLists.size() || it == _end)
    {
      return end();
    }

    return G4ManyFastLists_iterator<OBJECT>(trackList_it,
//                                           fAssociatedLists.begin(),
                                            it,
                                           &fAssociatedLists);
  }

template<class OBJECT>
typename G4ManyFastLists<OBJECT>::iterator G4ManyFastLists<OBJECT>::end()
  {
    if (fAssociatedLists.empty())
    {
      return G4ManyFastLists_iterator<OBJECT>(G4FastList_iterator<OBJECT>(),
                                             fAssociatedLists.end(),
                                             &fAssociatedLists);
    }

    return G4ManyFastLists_iterator<OBJECT>((fAssociatedLists.end()--)->end(),
                                           fAssociatedLists.end(),
                                           &fAssociatedLists);
  }

#include "G4ManyFastLists.icc"
#endif /* G4MANYFASTLISTS_HH_ */
