// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HadronicInteraction.cc,v 1.2 1999-12-15 14:52:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Interaction  base class
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 21-Mar-1997
 // Last modified: 04-Apr-1997
 
#include "G4HadronicInteraction.hh"
 
 G4double
  G4HadronicInteraction::GetMinEnergy(
   const G4Material *aMaterial, const G4Element *anElement ) const
  {
    G4int i;
    if( IsBlocked(aMaterial) )return 0.*GeV;
    if( IsBlocked(anElement) )return 0.*GeV;
    for( i=0; i<theMinCounterElements; ++i )
    {
      if( anElement == theMinElements[i] )return theMinEnergyListElements[i];
    }
    for( i=0; i<theMinCounter; ++i )
    {
      if( aMaterial == theMinMaterials[i] )return theMinEnergyList[i];
    }
    if( verboseLevel > 0 )
      G4cout << "*** Warning from HadronicInteraction::GetMinEnergy" << G4endl
           << "    material " << aMaterial->GetName()
           << " not found in min energy List" << G4endl;
    return theMinEnergy;
  }
 
 void
  G4HadronicInteraction::SetMinEnergy(
   G4double anEnergy,
   G4Element *anElement )
  {
    if( IsBlocked(anElement) )
      G4cout << "*** Warning from HadronicInteraction::SetMinEnergy" << G4endl
           << "    The model is not active for the Element  "
           << anElement->GetName() << "." << G4endl;
    
    for( G4int i=0; i<theMinCounterElements; ++i )
    {
      if( anElement == theMinElements[i] )
      {
        theMinEnergyListElements[i] = anEnergy;
        return;
      }
    }
    if( theMinCounterElements == MAX_LIST_SIZE )
      G4Exception("SetMinEnergy: exceeded size of min energy element List");
    theMinElements[theMinCounterElements] = anElement;
    theMinEnergyListElements[theMinCounterElements++] = anEnergy;
  }
 
 void
  G4HadronicInteraction::SetMinEnergy(
   G4double anEnergy,
   G4Material *aMaterial )
  {
    if( IsBlocked(aMaterial) )
      G4cout << "*** Warning from HadronicInteraction::SetMinEnergy" << G4endl
           << "    The model is not active for the Material "
           << aMaterial->GetName() << "." << G4endl;
    
    for( G4int i=0; i<theMinCounter; ++i )
    {
      if( aMaterial == theMinMaterials[i] )
      {
        theMinEnergyList[i] = anEnergy;
        return;
      }
    }
    if( theMinCounter == MAX_LIST_SIZE )
      G4Exception("SetMinEnergy: exceeded size of min energy material List");
    theMinMaterials[theMinCounter] = aMaterial;
    theMinEnergyList[theMinCounter++] = anEnergy;
  }
 
 G4double
  G4HadronicInteraction::GetMaxEnergy(
   const G4Material *aMaterial, const G4Element *anElement ) const
  {
    G4int i;
    if( IsBlocked(aMaterial) )return 0.0*GeV;
    if( IsBlocked(anElement) )return 0.0*GeV;
    for( i=0; i<theMaxCounterElements; ++i )
    {
      if( anElement == theMaxElements[i] )return theMaxEnergyListElements[i];
    }
    for( i=0; i<theMaxCounter; ++i )
    {
      if( aMaterial == theMaxMaterials[i] )return theMaxEnergyList[i];
    }
    if( verboseLevel > 0 )
      G4cout << "*** Warning from HadronicInteraction::GetMaxEnergy" << G4endl
           << "    material " << aMaterial->GetName()
           << " not found in min energy List" << G4endl;
    
    return theMaxEnergy;
  }
 
 void
  G4HadronicInteraction::SetMaxEnergy(
   G4double anEnergy,
   G4Element *anElement ) 
  {
    if( IsBlocked(anElement) )
      G4cout << "*** Warning from HadronicInteraction::SetMaxEnergy" << G4endl
           << "Warning: The model is not active for the Element  "
           << anElement->GetName() << "." << G4endl;
    
    for( G4int i=0; i<theMaxCounterElements; ++i )
    {
      if( anElement == theMaxElements[i] )
      {
        theMaxEnergyListElements[i] = anEnergy;
        return;
      }
    }
    if( theMaxCounterElements == MAX_LIST_SIZE )
      G4Exception("SetMaxEnergy: exceeded size of max energy element List");
    theMaxElements[theMaxCounterElements] = anElement;
    theMaxEnergyListElements[theMaxCounterElements++] = anEnergy;
  }

 void
  G4HadronicInteraction::SetMaxEnergy(
   G4double anEnergy,
   G4Material *aMaterial )
  {
    if( IsBlocked(aMaterial) ) 
      G4cout << "*** Warning from HadronicInteraction::SetMaxEnergy" << G4endl
           << "Warning: The model is not active for the Material "
           << aMaterial->GetName() << "." << G4endl;
    
    for( G4int i=0; i<theMaxCounter; ++i )
    {
      if( aMaterial == theMaxMaterials[i] )
      {
        theMaxEnergyList[i] = anEnergy;
        return;
      }
    }
    if( theMaxCounter == MAX_LIST_SIZE )
      G4Exception("SetMaxEnergy: exceeded size of max energy material List");
    theMaxMaterials[theMaxCounter] = aMaterial;
    theMaxEnergyList[theMaxCounter++] = anEnergy;
  }

 void 
  G4HadronicInteraction::DeActivateFor( G4Material *aMaterial )
  {
    if( theBlockedCounter == MAX_LIST_SIZE )
      G4Exception("DeActivateFor: exceeded size of blocked material List");
    theBlockedList[ theBlockedCounter++ ] = aMaterial;
  }

 void 
  G4HadronicInteraction::DeActivateFor( G4Element *anElement )
  {
    if( theBlockedCounterElements == MAX_LIST_SIZE )
      G4Exception("DeActivateFor: exceeded size of blocked elements List");
    theBlockedListElements[ theBlockedCounterElements++ ] = anElement;
  }

 G4bool 
  G4HadronicInteraction::IsBlocked( const G4Material *aMaterial ) const
  {
    G4bool tt = false;
    for( G4int i=0; i<theBlockedCounter; ++i )
    {
      if( aMaterial == theBlockedList[i] )
      {
        tt = true;
        break;
      }
    }
    return tt;
  }
 
 G4bool 
  G4HadronicInteraction::IsBlocked( const G4Element *anElement ) const
  {
    G4bool tt = false;
    for( G4int i=0; i<theBlockedCounterElements; ++i )
    {
      if( anElement == theBlockedListElements[i] )
      {
        tt = true;
        break;
      }
    }
    return tt;
  }
 
 /* end of file */
 
