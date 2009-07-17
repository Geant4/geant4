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
// $Id: G4InteractionContent.hh,v 1.5 2009-07-17 12:36:41 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4InteractionContent_h
#define G4InteractionContent_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4InteractionContent----------------
//             by Gunter Folger, June 1998.
//       class for a storing colliding particles in PartonString Models
// ------------------------------------------------------------

#include "globals.hh"

#include "G4VSplitableHadron.hh"
#include "G4Nucleon.hh"                // Uzhi 16.07.09
class G4InteractionContent 
{

  public:

      G4InteractionContent(G4VSplitableHadron *aPrimaryParticipant);

      ~G4InteractionContent();

      int operator==(const G4InteractionContent &right) const;
      int operator!=(const G4InteractionContent &right) const;
      
      G4VSplitableHadron * GetProjectile() const ;
      G4VSplitableHadron * GetTarget() const;

      void                 SetTargetNucleon(G4Nucleon * aNucleon); // Uzhi 16.07.09
      G4Nucleon          * GetTargetNucleon() const;              // Uzhi 16.07.09

      void SetTarget(G4VSplitableHadron *aTarget);

      G4int GetNumberOfSoftCollisions();
      G4int GetNumberOfHardCollisions();
      void  SetNumberOfSoftCollisions(int);
      void  SetNumberOfHardCollisions(int);
      G4int GetNumberOfDiffractiveCollisions();
      void  SetNumberOfDiffractiveCollisions(int);

      void SplitHadrons();
 
#ifdef G4DEBUG
      void Dump();
#endif      
     

public:
      G4InteractionContent(){}
      G4InteractionContent(const G4InteractionContent &right);
      const G4InteractionContent & operator=(const G4InteractionContent &right);

  protected:

  private:

      G4VSplitableHadron * theTarget;
      G4VSplitableHadron * theProjectile;
      G4Nucleon          * theTargetNucleon;
      
      G4int theNumberOfHard;
      G4int theNumberOfSoft;
      G4int theNumberOfDiffractive;
};

// Class G4InteractionContent 

inline G4VSplitableHadron * G4InteractionContent::GetProjectile() const
{
	return theProjectile;
}

inline G4VSplitableHadron * G4InteractionContent::GetTarget() const
{
	return theTarget;
}

inline void G4InteractionContent::SetTarget(G4VSplitableHadron *aTarget)
{
	theTarget = aTarget;
}

inline void G4InteractionContent::SetTargetNucleon(G4Nucleon * aNucleon) // Uzhi 16.07.09
{
        theTargetNucleon = aNucleon;
}

inline G4Nucleon * G4InteractionContent::GetTargetNucleon() const       // Uzhi 16.07.09
{
       return theTargetNucleon;
}

inline G4int G4InteractionContent::GetNumberOfSoftCollisions()
{
	return theNumberOfSoft;
}

inline G4int G4InteractionContent::GetNumberOfHardCollisions()
{
	return theNumberOfHard;
}

inline void G4InteractionContent::SetNumberOfSoftCollisions(int nCol)
{
	theNumberOfSoft = nCol;
}

inline void G4InteractionContent::SetNumberOfHardCollisions(int nCol)
{
	theNumberOfHard = nCol;
}

inline G4int G4InteractionContent::GetNumberOfDiffractiveCollisions()
{
	return theNumberOfDiffractive;
}

inline void G4InteractionContent::SetNumberOfDiffractiveCollisions(int nCol)
{
	theNumberOfDiffractive = nCol;
}

inline void G4InteractionContent::SplitHadrons()
{
	if ( theProjectile != NULL ) theProjectile->SplitUp();
	if ( theTarget != NULL ) theTarget->SplitUp();
}
#ifdef G4DEBUG
inline void G4InteractionContent::Dump()
{
	G4cout  << " G4InteractionContent " << this << G4endl
	        << "Hard/Soft/Diff " 
		<< theNumberOfHard<<" / " 
		<<theNumberOfSoft<<" / " 
		<<theNumberOfDiffractive << G4endl
		<< "Projectile " ;
        if ( theProjectile ) 
	{ G4cout <<  theProjectile->GetDefinition()->GetPDGEncoding() 
		 << "  " << theProjectile->Get4Momentum()<< G4endl;
	} else {	 
	  G4cout << " none " << G4endl;
	}	    
        if ( theTarget ) 
	{ G4cout <<  theTarget->GetDefinition()->GetPDGEncoding() 
		 << "  " << theTarget->Get4Momentum()<< G4endl;
	} else {	 
	  G4cout << " none " << G4endl;
	}	    
		    
}      
#endif      

#endif


