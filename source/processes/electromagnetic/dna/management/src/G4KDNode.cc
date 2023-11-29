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
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4KDNode.hh"
#include "G4KDTree.hh"
#include <ostream>

//*********************************************

//______________________________________________________________________
// Node functions

//void* GetData(G4KDNode* node)
//{
//	return node->GetData() ;
//}
//
//const double* GetNodePosition(G4KDNode* node)
//{
//	return node->GetPosition() ;
//}

//______________________________________________________________________

void InactiveNode(G4KDNode_Base* node)
{
  if(node == nullptr) return;
//	if(node->IsValid())
  node->InactiveNode();
}

void Free(G4KDNode_Base*& node)
{
  if(node) delete node;
  node = nullptr;
}

//______________________________________________________________________
G4KDNode_Base::G4KDNode_Base(G4KDTree* tree,
		G4KDNode_Base* parent):
    						fTree(tree),
    						fLeft(0), fRight(0), fParent(parent)
{
  fSide = 0;
  fAxis = fParent == 0? 0 : fParent->fAxis +1 < fTree->fDim? fParent->fAxis+1:0;
}

// Copy constructor should not be used
G4KDNode_Base::G4KDNode_Base(const G4KDNode_Base& ):
    						fTree(0),
    						fLeft(0), fRight(0), fParent(0)
{
  fSide = 0;
  fAxis = 0;
}

// Assignement should not be used
G4KDNode_Base& G4KDNode_Base::operator=(const G4KDNode_Base& right)
{
	if (this == &right) return *this;
	fTree = right.fTree;
	fLeft = right.fLeft;
	fRight = right.fRight;
	fParent = right.fParent;
	fSide = right.fSide;
	fAxis = right.fAxis;
	return *this;
}

G4KDNode_Base::~G4KDNode_Base()
{
}

void G4KDNode_Base::InactiveNode()
{
	fTree->NoticeNodeDeactivation();
}

G4int G4KDNode_Base::GetDim() const
{
	if(fTree)
          return (G4int)fTree->GetDim();
	else
          return -1;
}

G4int G4KDNode_Base::Insert(G4KDNode_Base* newNode)
{
	G4KDNode_Base* aParent = FindParent(*newNode);
	// TODO check p == aParent->pos
	// Exception

	newNode->fAxis = aParent->fAxis +1 < fTree->GetDim()? aParent->fAxis+1:0;
	newNode->fParent = aParent ;

	if((*newNode)[aParent->fAxis] > (*aParent)[aParent->fAxis])
	{
		aParent->fRight = newNode ;
		newNode->fSide = 1 ;
	}
	else
	{
		aParent->fLeft = newNode ;
		newNode->fSide = -1 ;
	}

	newNode->fRight = 0;
	newNode->fLeft = 0;

	return 0 ;
}


void G4KDNode_Base::PullSubTree()
{
	if(fParent)
	{
		if(fSide == -1)
		{
			fParent->fLeft = 0;
		}
		else
			fParent->fRight = 0;
	}
	if(fLeft) fLeft -> PullSubTree();
	if(fRight) fRight-> PullSubTree();

	fParent  = 0 ;
	fRight   = 0 ;
	fLeft    = 0 ;
	fTree    = 0 ;
}

void G4KDNode_Base::RetrieveNodeList(std::list<G4KDNode_Base *>& output)
{
	output.push_back(this);

	if(fLeft)
		fLeft->RetrieveNodeList(output);

	if(fRight)
		fRight->RetrieveNodeList(output);
}

void G4KDNode_Base::Print(std::ostream& out, int level) const
{
	// Print node level
	out << G4endl;
	for (G4int i=0; i<level; ++i)         // Indent to level
	{
		out << "  ";
	}
	out << level;

	// Print children
	if(fLeft)
	{
		fLeft->Print(out, level + 1);
	}
	if(fRight)
	{
		fRight->Print(out, level + 1);
	}
}
