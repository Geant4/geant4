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
// $Id: tree.hh,v 1.2 2002-06-04 09:23:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Tree
//
// Very simple tree data-type to get rid of QT-ListView ...
// Trees are build once and not to be modified.
// No memory management for Data!
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapTree_h
#define OlapTree_h

#include "g4std/set"
#include "g4std/vector"
#include "g4std/list"

#include "G4Types.hh"

template <class Data> class TreeNodeIterator;

template <class Data>
class TreeNode
{

  friend class TreeNodeIterator<Data>;
  
public:

  TreeNode(TreeNode * parent, Data data)
    : parent_(parent), firstChild_(0),
      lastChild_(0), nextSibling_(0), data_(data)
    { 
      if (parent)
      {
        if(!parent->firstChild_)
        {
	  parent->firstChild_=this;
	  parent->lastChild_=this;
	}  
	else
        {
	  parent->lastChild_->nextSibling_ = this;
	  parent->lastChild_=this;
	}   
      }	
    }
    
  ~TreeNode() // deletes the whole subtree as well!
    { 
      if (parent_)
      {
        TreeNode * i = parent_->firstChild_;
	while(*i)
        {
	  TreeNode * temp = i;
	  i = i->nextSibling_;
	  delete temp;
	}
      } 
    } 

  TreeNode * parent() const { return parent_; }
  Data data() const { return data_; }
  
  G4int childCount() const 
  {
    TreeNode * i = firstChild_;
    G4int c=0;
    while (i)
    {
      i = i->nextSibling_; c++;
    }
    return c;  
  }
  TreeNode * firstChild() const { return firstChild_; }
  TreeNode * nextSibling() const { return nextSibling_; }
  TreeNode * lastChild() const { return lastChild_; }  
  
protected:

  TreeNode * parent_;
  TreeNode * firstChild_;
  TreeNode * lastChild_;
  TreeNode * nextSibling_;
  Data data_;
  
private:  
   TreeNode(const TreeNode &);
   TreeNode & operator=(const TreeNode *);
};


template <class Data>
class Tree
{
public:

  typedef TreeNode<Data> node_t;

  Tree(Data rt) 
    : root_(new node_t(0,rt)) {}
    
  node_t * root() { return root_; }
  
protected:

  node_t * root_;
};


template<class Data>
class TreeNodeIterator
{
public:
  TreeNodeIterator(TreeNode<Data> * np)
   : myroot_(np), curpos_(np)
     { 
       //if (np->parent())
         //nexts_=++(np->parent()->firstChild());
     }
   
  ~TreeNodeIterator() { }
  
  TreeNode<Data> * current() { return curpos_; }

  TreeNode<Data> * next()
  {
     
     if (curpos_->firstChild_)
     {
       curpos_ = curpos_->firstChild_;
       return curpos_;
     }  
     
     if (curpos_->nextSibling_)
     {
       curpos_ = curpos_->nextSibling_;
       return curpos_;
     }  
     while(up())
     {
       if (curpos_->nextSibling_)
       {
         curpos_ = curpos_->nextSibling_;
         return curpos_;	 
       }	 
     }  
     return 0;   	     
  }
  
protected:

  G4bool up()
  {
    G4bool result(true);
    if (curpos_->parent_)
      curpos_ = curpos_->parent_;
    else
      result=false;
    
    if (curpos_==myroot_)
    {
      result=false;  
    }   
    return result;    
  } 
  
  TreeNode<Data> * myroot_;  
  TreeNode<Data> * curpos_;
  //Tree<Data>::c_iterator nextc_;
};

#endif
