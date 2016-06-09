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
// $Id: G4KDTreeResult.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4KDTREERESULT_HH
#define G4KDTREERESULT_HH

#include <list>
#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"
class G4KDTree;
class G4KDNode;
struct ResNode;

class G4KDTreeResult;

typedef G4ReferenceCountedHandle<G4KDTreeResult> G4KDTreeResultHandle;
typedef G4ReferenceCountedHandle<ResNode> ResNodeHandle;

/**
  * G4KDTreeResult enables to go through the nearest entities found
  * by G4KDTree.
  */

class G4KDTreeResult : protected std::list<ResNode>
{
protected :
    G4KDTree *fTree;
    std::list<ResNode>::iterator fIterator;

public:
    G4KDTreeResult(G4KDTree*);
    virtual ~G4KDTreeResult();

    void Insert(double, G4KDNode*);

    void Clear();

    void Sort();

    /* returns the size of the result set (in elements) */
    size_t GetSize();

    size_t size();

    /* rewinds the result set iterator */
    void Rewind();

    /* returns non-zero if the set iterator reached the end after the last element */
    bool End();

    /* advances the result set iterator
     */
    void Next();

    /* returns the data pointer (can be null) of the current result set item
     * and optionally sets its position to the pointers(s) if not null.
     */
    void* GetItemData();
    void* GetItem(double*& /*position*/);
    void* GetItem(double& x, double& y, double& z); // 3D
    void* GetItemNDistanceSQ(double& /*distance*/);
    void* GetItemNDistanceSQ(double*& /*position*/, double& /*distance*/);
    double GetDistanceSqr();
};

#endif // G4KDTREERESULT_HH
