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
/// \file geometry/olap/src/G4GeoNav.cc
/// \brief Implementation of the G4GeoNav class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// G4GeoNav
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "G4GeoNav.hh"

G4GeoNav::G4GeoNav( G4LogicalVolume * rv )
  : theLVTree(new LVTree(rv))
{
   theCurLV = theLVTree->root();
   theRootLV = theLVTree->root();
   RecursiveFill(theLVTree->root());
}  

G4GeoNav::~G4GeoNav()
{}


G4LogicalVolume * G4GeoNav::NextLV()
{
  //FIXME: G4GeoNav::NextLV() 

  TreeNodeIterator<G4LogicalVolume*> tni(theCurLV); 

  G4LogicalVolume * v = 0;
  if (tni.next()) {
    v = tni.current()->data();
    theCurLV = tni.current();
  }  
  return v;  
/*        
  G4cerr << "DON'T CALL G4GeoNav::Next()!" << endl;
  exit(1);
  return 0;
*/
}


G4int G4GeoNav::FilterLV(const G4String & aRegexStr,
                         std::vector<G4LogicalVolume*> & result, 
                         G4bool stopAtFirst)
{
  regex_t aRegex;
  const char * aRegexCStr = aRegexStr.data();
  if (regcomp(&aRegex,aRegexCStr,0)) {
    G4cerr << "failed to interpret regex-string" << G4endl;
    return 0;
  }
  
  LVTree::node_t * aNode = theLVTree->root();
  FindLV(&aRegex,aNode,result,stopAtFirst);
  regfree(&aRegex);
  return result.size();  
}                         


void G4GeoNav::FindLV(regex_t * aRegex, LVTree::node_t * node, 
                      std::vector<G4LogicalVolume*>& result,
                      G4bool stopAtFirst)
{
  if( !regexec(aRegex, node->data()->GetName().data(), 0,0,0))
  { 
    result.push_back(node->data());
    if (stopAtFirst)
    {
      theCurLV = node;
      return;
    } 
  }
  
  LVTree::node_t * i = node->firstChild();
  while(i) { // recursive
    FindLV(aRegex, i, result, stopAtFirst);
    i = i->nextSibling();
  }  

}        

G4int G4GeoNav::PathLV(std::vector<G4LogicalVolume*> & result)
{
   result.push_back(theCurLV->data());
   G4int level=1;
   LVTree::node_t * node = theCurLV;
   while (node->parent()) 
   {
     node = node->parent();
     G4LogicalVolume * aLV = node->data();
     result.push_back(aLV);
     level++;
   }
   return level;  
}
                         
                      
G4int G4GeoNav::Tokenize(const G4String & aStr,
                         std::vector<G4String>& tokens)
{
    G4String::size_type c = aStr.size();
    G4String::size_type idx = 0;
    G4String::size_type idx2 = 0;
    
    // empty string means 'root'
    if (!c)
    {
       tokens.push_back("/");
       return 1;
    }   
          
    // scan for '/'
    std::vector<G4int> pos;
    pos.push_back(0); // begin of input  
    G4String sep("/");
       
    if (aStr[idx]==sep[idx]) 
      tokens.push_back("/");
     
    G4String curStr = aStr; // G4String("-")  ;  

    G4String::size_type i=0;
    for (i=0; i<c; i++) {
       idx = i;
       if(curStr[idx]==sep[idx2]) {
          #ifdef QT_DEBUG_3
            G4cout << i << " " << curStr[i] << G4endl;
          #endif 
          pos.push_back(i);
       }  
    }    
    pos.push_back(c-1); // end of input
    
    #ifdef QT_DEBUG_3
       for (G4int i=0; i<pos.size(); i++ )
         G4cout << pos[i] << " ";
       G4cout << G4endl;
    #endif
     
    for (i=1; i<pos.size(); i++)
    {
       G4String newString;
                
       for(G4int j=pos[i-1]; j<=pos[i]; j++) {
          idx = j;
          if (curStr[idx]!=sep[idx2])
             newString.append(curStr[idx]);
       }
       if (newString.size())
         tokens.push_back(newString);
    }    
    
    #ifdef QT_DEBUG 
    for (G4int i=0; i<tokens.size(); i++) 
      G4cout << "tok " << i << " " << tokens[i] << G4endl;
    #endif                                          
    
    return tokens.size();
}


G4LogicalVolume * G4GeoNav::ChangeLV(const G4String & aRegExp)
{
    std::vector<G4String> tokens;
    G4int c = Tokenize(aRegExp,tokens);
    
    LVTree::node_t * anItem = theCurLV;
    LVTree::node_t * anItemBefore = theCurLV;   
    if (c) {
      
      std::vector<G4String>::iterator it = tokens.begin();
      if (*it==G4String("/"))
      {                        // absolute or relative
         anItem = theRootLV;
         it++;
      } 
      
      while ( it != tokens.end() )
      {
        //regex_t aRegex;
        if (*it == G4String(".."))
        {
           anItem = anItem->parent();
           if (!anItem) {
              G4cerr << "not found!" << G4endl;
              return 0;
           }     
           it++;    
        }
        else
        {
           regex_t aRegex;
           const char * aRegexCStr = (*it).data();
           if (regcomp(&aRegex,aRegexCStr,0))
           { 
               G4cerr << "failed to interpret regex-string" << G4endl;
               return 0;
           }           
           
           G4int children = anItem->childCount();
           
           if (!children)
           {
              G4cerr << "not found!!" << G4endl;
              return 0;
           }
           
           //LVTree::node_t * temp = *(anItem->firstChild());
           //anItem = temp;
           
           // loop over children
           LVTree::node_t * u= anItem->firstChild();
           G4bool found=false;
           while (u)
           {   
             if ( !regexec(&aRegex, u->data()->GetName().data(),0,0,0))
             {
                found=true;
                anItem = u;
                break;   
             }
             u = u->nextSibling();
           }
           
           regfree(&aRegex);
           if (!found)
           {
              G4cerr << "not found!!!!" << G4endl;
              return 0;
           }   
           it++;        
        }
        
      }
      
      if (anItem!=anItemBefore)
      {
        theCurLV=anItem;
        return anItem->data();
      }         
      return 0;            
    }
  return 0;  
}

G4int G4GeoNav::PwdLV(std::vector<G4LogicalVolume *>& result)
{
   result.push_back(theCurLV->data());
   LVTree::node_t * temp = theCurLV;
   while (temp->parent())
   {
     temp = temp->parent();
     result.push_back(temp->data());
   }
   return result.size();  
}


void G4GeoNav::RecursiveFill(LVTree::node_t* node)
{
   G4LogicalVolume * lv = node->data();
   
   std::set<G4LogicalVolume*> lvset;
   
   for ( G4int i=0; i<lv->GetNoDaughters(); ++i)
   {
     lvset.insert(lv->GetDaughter(i)->GetLogicalVolume());  
   }
   
   std::set<G4LogicalVolume*>::iterator it = lvset.begin();
   for(;it!=lvset.end();++it)
   {
      LVTree::node_t * newnode = new LVTree::node_t(node, *it);
        //create new node '*it' in parent 'node'
      RecursiveFill(newnode);
   }  

}

G4int G4GeoNav::LsLV(std::vector<G4LogicalVolume *>& result)
{
  G4int c=0;
  LVTree::node_t * n = theCurLV->firstChild();
  while(n)
  {
    result.push_back(n->data());
    n = n->nextSibling();
    c++;
  }   
  return c;
}
