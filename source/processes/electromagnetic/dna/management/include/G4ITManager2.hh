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
// $Id: G4ITManager.hh 79462 2014-03-01 17:15:17Z matkara $
//
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source, article citations are crucial.
// If you use Geant4-DNA chemistry and you publish papers about your software, in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 
//----------------------------------------------------------------


#ifdef DNADEV
#ifndef G4ITManager2_hh
#define G4ITManager2_hh 1

#include "globals.hh"
#include <map>
#include "G4AllITManager.hh"
#include "G4ITBox.hh"
#include "G4KDTree.hh"
#include "G4Track.hh"
#include "G4TrackList.hh"
#include "G4memory.hh"

// namespace
enum EListType
{
	ITMain = 0,
	ITSecondaries = 1,
	ITWaiting= 2,
	ITDelayed= 3,
	ITKill= 4,
	ITUndefined= -1
};

struct Lists
{
	G4TrackList* GetList(EListType type)
	{
		switch(type)
		{
		case ITMain:
		{
			return fpMainList;
			break;
		}
		case ITSecondaries:
		{
			return &fSecondaries;
			break;
		}
		case ITWaiting:
		{
			return fpWaitingList; break;
		}
		case ITKill:
		{
			return &fToBeKilledList; break;
		}
		case ITDelayed:
		{
			return 0;
		}
		default:
		{
			return 0; break;
		}
		}
		return 0;
	}

	G4TrackList* GetMainList()
	{
		return fpMainList;
	}

	G4TrackList* GetListOfSecondaries()
	{
		return &fSecondaries;
	}

	G4TrackList* GetWaitingList()
	{
		return fpWaitingList;
	}

	G4TrackList* GetListToBeKilled()
	{
		return &fToBeKilledList;
	}

	std::map<double,G4TrackList* >& GetDelayedList(){
		return fDelayedList;
	}

	Lists()
	{
		fpMainList = new G4TrackList();
		fpWaitingList = 0;
	}

	Lists(Lists& right)
	{
		fpMainList = right.fpMainList;
		fpWaitingList = right.fpWaitingList;
		right.fSecondaries.transferTo(&fSecondaries);
		right.fDelayedList.swap(fDelayedList);
	}

	~Lists()
	{
		if(fpMainList)
		{
			delete fpMainList;
		}

		if(fpWaitingList)
		{
			delete fpWaitingList;
		}
	}

	void push(G4Track* track, EListType type)
	{
		if(type == ITDelayed)
		{
			pushToDelayed(track);
		}
		else
		{
			GetList(type)->push_back(track);
		}
	}

	void pop(G4Track* track, EListType type)
	{
		if(type != ITDelayed)
		{
			GetList(type)->pop(track);
		}
	}

	void pushToDelayed(G4Track* track)
	{
		double time = track->GetGlobalTime();
		std::map<double,G4TrackList* >::iterator it = fDelayedList.find(time);

		if(it == fDelayedList.end())
		{
			std::pair<std::map<double,G4TrackList* >::iterator,bool> ret =
					fDelayedList.insert(std::pair<double,G4TrackList* >(time,new G4TrackList()));
			ret.first->second->push_back(track);
		}
		else
		{
			it->second->push_back(track);
		}
	}

	G4TrackList*                    fpMainList ;
	G4TrackList                     fSecondaries ; //to merge with fpMainList
	G4TrackList*                    fpWaitingList ; // Waiting queue of currentList
	std::map<double,G4TrackList* >  fDelayedList ;
	G4TrackList                     fToBeKilledList ;
};

class ListHandler
{
public :
	G4shared_ptr<Lists> fLists;

	ListHandler()
	{}

	Lists* operator->()
	{
		if(!fLists)
		{
			fLists.reset(new Lists());
		}

		return fLists.get();
	}

	Lists* operator*()
	{
		if(!fLists)
		{
			fLists.reset( new Lists() );
		}

		return fLists.get();
	}
};

class ListGroup
{
public:

	typedef std::map<G4ITType, std::map<G4ITType, ListHandler> > MapMapList;
	typedef std::map<G4ITType, ListHandler> MapList;

public:
	MapMapList fITList ;
	MapList* fLastListMap;

public:
	class iterator
	{
	public :
		MapList* fMapLists;
		MapList::iterator fMapIterator;
		EListType fListType;
		G4TrackList* fpCurrentTrackList;
		G4TrackList::iterator fTrackIterator;

		iterator(MapList* list, EListType type, bool end = false)
		{
			fMapLists = list;
			fListType = type;
			if(end == false)
			{
				fMapIterator = fMapLists->begin();
				fpCurrentTrackList = GetCurrentTrackList();
				if(fpCurrentTrackList)
				{
					fTrackIterator = fpCurrentTrackList->begin();
				}
			}
			else
			{
				fpCurrentTrackList = 0;
				fMapIterator = fMapLists->end();
			}
		}

		iterator() : fMapLists(0)
		{
			fListType = ITUndefined;
			fpCurrentTrackList = 0;
		}

		iterator(const iterator& right)
		{
			fMapLists = right.fMapLists;
			fMapIterator = right.fMapIterator;
			fListType = right.fListType;
			fpCurrentTrackList = right.fpCurrentTrackList ;
			fTrackIterator = right.fTrackIterator;
		}

		iterator& operator=(const iterator& right)
		{
			if(this == &right) return *this;
			fMapLists = right.fMapLists;
			fMapIterator = right.fMapIterator;
			fListType = right.fListType;
			fpCurrentTrackList = right.fpCurrentTrackList ;
			fTrackIterator = right.fTrackIterator;
			return *this;
		}

		bool operator==(const iterator& right) const{

			if(fpCurrentTrackList == 0 ) //&& right.fpCurrentTrackList == 0)
			{

				G4cout << "both fpCurrentTrackList = 0 " << G4endl;
				return (
						fMapLists == right.fMapLists &&
						fMapIterator == right.fMapIterator &&
						fListType == right.fListType);
			}

			return (
					fMapLists == right.fMapLists &&
					fMapIterator == right.fMapIterator &&
					fListType == right.fListType &&
					fTrackIterator == right.fTrackIterator);
		}

		bool operator!=(const iterator& right) const{
			return !(*this == right);
		}

		G4Track* operator->()
		{
			return *fTrackIterator;
		}

		const G4Track* operator->() const
		{
			return *fTrackIterator;
		}

		G4Track* operator*()
		{
			return *fTrackIterator;
		}

		const G4Track* operator*() const
		{
			return *fTrackIterator;
		}

		G4TrackList* GetCurrentTrackList()
		{
			if(fMapIterator == fMapLists->end())
			{
				G4cout << "end !! " << G4endl;
				return 0;
			}
			return fMapIterator->second->GetList(fListType);
		}

		iterator& operator++(int) {
			if(fMapIterator == fMapLists->end()) return *this;

			fTrackIterator++;
			if(fTrackIterator == fpCurrentTrackList->end())
			{
				while(
						GetCurrentTrackList() &&
						fTrackIterator == fpCurrentTrackList->end()
				)
				{
					fMapIterator++;
					fpCurrentTrackList = GetCurrentTrackList();
					if(fpCurrentTrackList)
					{
						fTrackIterator = fpCurrentTrackList->begin();
					}
					else
					{
						G4cout << "break" << G4endl;
						break;
					}
				};

				return *this;
			}
			else
			{
				return *this;
			}
			return *this;
		}

		bool operator()()
		{
			(*this)++;
			return (fMapIterator != fMapLists->end());
		}

	};

	class all_iterator
	{
		ListGroup* fGroup;
		MapMapList* fMapMapLists;
		G4ITType fType;
		EListType fListType;
		MapMapList::iterator fMapMapListIterator;
		iterator ListIterator;

	public:
		all_iterator(ListGroup* group,
				G4ITType ITType,
				EListType listType,
				bool end = false) {
			fGroup = group;
			fMapMapLists = &(group->fITList);
			fType = ITType;
			fListType = listType;
			if(end == false)
			{
				fMapMapListIterator = fMapMapLists->begin();
				ListIterator = iterator(&(fMapMapLists->begin()->second), listType, end);
			}
			else
			{
				fMapMapListIterator = fMapMapLists->end();
				ListIterator = iterator(&(fMapMapLists->begin()->second), listType, end);
			}
		}

		all_iterator(const all_iterator& right) {
			fGroup = right.fGroup;
			fMapMapLists = right.fMapMapLists;
			fType = right.fType;
			fListType = right.fListType;
			fMapMapListIterator = right.fMapMapListIterator;
			ListIterator = right.ListIterator;
		}

		bool operator==(const all_iterator& right) {

			if(fGroup == right.fGroup &&
					fMapMapLists  ==  right.fMapMapLists &&
					fType == right.fType &&
					fListType == right.fListType
			)
			{
				if(fMapMapListIterator == fMapMapLists->end() && right.fMapMapListIterator == fMapMapLists->end())
				{
					G4cout << "Return true" << G4endl;
					return true;
				}

				G4cout << "Return " <<  (fMapMapListIterator ==  right.fMapMapListIterator &&
						ListIterator == right.ListIterator) << G4endl;

				return (
						fMapMapListIterator ==  right.fMapMapListIterator &&
						ListIterator == right.ListIterator
				);

			}

			G4cout << "Return false" << G4endl;
			return false;
		}

		bool operator!=(const all_iterator& right) {
			return !(*this == right);
		}

		G4Track* operator*()
		{
			return *ListIterator;
		}

		MapList* GetCurrentMapList()
		{
			if(fMapMapListIterator == fMapMapLists->end()) return 0;
			return &(fMapMapListIterator->second);
		}

		all_iterator& operator++(int) {

			G4cout << "all_iterator& operator++" << G4endl;

			if(GetCurrentMapList() == 0)
			{
				G4cout << "all_iterator& operator++ : break" << G4endl;
				return *this;
			}

			while(ListIterator() == false &&
					GetCurrentMapList())
			{
				G4cout << "in while " << G4endl;
				fMapMapListIterator++;
				if(GetCurrentMapList())
				{
					ListIterator = iterator(GetCurrentMapList(), fListType, false);
				}
				else
				{
					G4cout << "all_iterator : break " << G4endl;
					break;
				}
			}
			return *this;
		}
		/*
		bool operator()()
		{
			(*this)++;
		}
		*/
	};

public:
	void push(G4Track* track, EListType type)
	{
		G4IT* IT = GetIT(track);
		fITList[IT->GetITType()][IT->GetITSubType()]->push(IT->GetTrack(),type);
	}

	void pop(G4Track* track, EListType type)
	{
		G4IT* IT = GetIT(track);
		MapMapList::iterator it = fITList.find(IT->GetITType());

		if(it == fITList.end())
		{
			return;
		}

		MapList::iterator it2 = it->second.find(IT->GetITSubType());

		if(it2 == it->second.end())
		{
			return;
		}

		it2->second->pop(track,type);
	}

	void Clear()
	{
		fITList.clear();
		fLastListMap = 0;
	}

	all_iterator begin(EListType type)
	{
		all_iterator it = all_iterator(this, -1, type, false);
		return it;
	}

	all_iterator end(EListType type)
	{
		all_iterator it = all_iterator(this, -1, type, true);
		return it;
	}

	iterator begin(G4ITType ITType, EListType type)
	{
		iterator it = iterator(&fITList[ITType], type);
		return it;
	}

	iterator end(G4ITType ITType, EListType type)
	{
		iterator it = iterator(&fITList[ITType], type, true);
		return it;
	}

	bool FindSubListMap(const G4ITType& ITType)
	{
		MapMapList::iterator it = fITList.find(ITType);
		if(it == fITList.end())
		{
			fLastListMap = 0;
			return 1;
		}
		fLastListMap = &(it->second);
		return 0;
	}

	G4TrackList* GetTrackList(const G4ITType& ITType, EListType type, G4ITType ITSubType = 0)
	{
		if(FindSubListMap(ITType)) return 0;

		MapList::iterator it2 = fLastListMap->find(ITSubType);
		if(it2 == fLastListMap->end())	return it2->second->GetList(type);

		return 0;
	}

	Lists* GetList(const G4ITType& ITType)
	{
		if(FindSubListMap(ITType)) return 0;

		std::map<G4ITType, ListHandler>::iterator it = fLastListMap->find(ITType);
		if(it != fLastListMap->end()) return *(it->second);

		return 0;
	}
};

/**
 * G4VITManager is just a virtual interface for G4ITManager.
 * For more details, please have a look at the description
 * of ITManager.
 */

class G4ITManager2
{
public :
	G4ITManager2()
{
		fVerbose= 0;
		gITManager2 = this;
}
	virtual ~G4ITManager2(){;}
	static G4ITManager2*   Instance()
	{
		if(gITManager2 == 0)
		{
			new G4ITManager2();
		}

		return gITManager2;
	}

	virtual void PushMain(G4Track* track)
	{
		fListGroup.push(track, ITMain);
	}

	virtual void Push(G4Track* track, EListType type)
	{
		fListGroup.push(track, type);
	}

	ListGroup* GetListGroup()
	{
		return &fListGroup;
	}

	virtual void PopMain(G4Track* track)
	{
		fListGroup.pop(track, ITMain);
	}

	virtual void Pop(G4Track* track, EListType type)
	{
		fListGroup.pop(track, type);
	}

	virtual void Clear()
	{
		fListGroup.Clear();
	}

	//Make a method to extract number of IT from Box
	G4int NbElements(const G4IT*);

	void SetVerboseLevel(G4int level)
	{
		fVerbose = level;
	}
	G4int GetVerboseLevel()
	{
		return fVerbose;
	}
	/*
	G4TrackList* GetTrackList(const G4IT* IT, EListType type, G4ITType subType = 0)
	{
		return fListGroup.GetTrackList(IT, type, subType);
	}
	*/
protected :
	G4ITType fType ;
	G4int fVerbose;
	ListGroup fListGroup;
	static G4ITManager2* gITManager2;
};

#endif
#endif
