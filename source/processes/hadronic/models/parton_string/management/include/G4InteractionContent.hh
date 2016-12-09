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
// $Id: G4InteractionContent.hh 100828 2016-11-02 15:25:59Z gcosmo $
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

//#include "G4PhysicalConstants.hh"       // Uzhi 8.09.2014
//#include "G4SystemOfUnits.hh"           // Uzhi 8.09.2014

#include "G4VSplitableHadron.hh"
#include "G4Nucleon.hh"

//#define debug_QGSM
#ifdef debug_QGSM
  #include "../../qgsm/include/G4QGSMSplitableHadron.hh"
  #include <deque>
#endif  

class G4InteractionContent 
{
  public:
    G4InteractionContent() {}
    G4InteractionContent(G4VSplitableHadron *aPrimaryParticipant);

    ~G4InteractionContent();

    G4bool operator<(const G4InteractionContent &right) const;
      
    G4VSplitableHadron * GetProjectile() const ;
    G4VSplitableHadron * GetTarget() const;

    void                 SetProjectileNucleon(G4Nucleon * aNucleon);
    G4Nucleon          * GetProjectileNucleon() const;

    void                 SetTargetNucleon(G4Nucleon * aNucleon);
    G4Nucleon          * GetTargetNucleon() const;

    void SetTarget(G4VSplitableHadron *aTarget);

    G4int GetNumberOfSoftCollisions();
    G4int GetNumberOfHardCollisions();
    void  SetNumberOfSoftCollisions(int);
    void  SetNumberOfHardCollisions(int);
    G4int GetNumberOfDiffractiveCollisions();
    void  SetNumberOfDiffractiveCollisions(int);

    void SplitHadrons();

    void     SetInteractionTime(G4double aValue);
    G4double GetInteractionTime() const;         
    void     SetStatus(G4int aValue);            
    G4int    GetStatus() const;                  
 
    #ifdef G4DEBUG
    void Dump();
    #endif

  private:
    G4InteractionContent & operator=(const G4InteractionContent &right);
    G4InteractionContent(const G4InteractionContent &right);
    int operator==(const G4InteractionContent &right) const;
    int operator!=(const G4InteractionContent &right) const;

  protected:

  private:
    G4VSplitableHadron * theTarget;
    G4VSplitableHadron * theProjectile;

    G4Nucleon          * theProjectileNucleon;
    G4Nucleon          * theTargetNucleon;
      
    G4int theNumberOfHard;
    G4int theNumberOfSoft;
    G4int theNumberOfDiffractive;

    G4double theInteractionTime;
    G4int    curStatus;
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

inline void G4InteractionContent::SetProjectileNucleon(G4Nucleon * aNucleon)
{
  theProjectileNucleon = aNucleon;
}

inline G4Nucleon * G4InteractionContent::GetProjectileNucleon() const
{
  return theProjectileNucleon;
}

inline void G4InteractionContent::SetTargetNucleon(G4Nucleon * aNucleon)
{
  theTargetNucleon = aNucleon;
}

inline G4Nucleon * G4InteractionContent::GetTargetNucleon() const
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
  //G4cout<<"InterContent Proj "<<theProjectile<<G4endl;
  if ( theProjectile != NULL ) {theProjectile->SplitUp();}
  //G4cout<<"InterContent Targ "<<theTarget<<G4endl;
  if ( theTarget != NULL )     {theTarget->SplitUp();}
  #ifdef G4DEBUG
  //Dump();
  #endif
}

#ifdef G4DEBUG
inline void G4InteractionContent::Dump()
{
  G4LorentzVector mom(0.,0.,0.,0.);
  G4cout  << " G4InteractionContent " << this << G4endl
	  << "Hard/Soft/Diff "
	  << theNumberOfHard<<" / "
	  <<theNumberOfSoft<<" / "
	  <<theNumberOfDiffractive << G4endl
	  << "Projectile " ;

  if ( theProjectile ) {
    G4cout <<  theProjectile->GetDefinition()->GetPDGEncoding()
           << "  " << theProjectile->Get4Momentum()<< G4endl;
    mom+=theProjectile->Get4Momentum();
    #ifdef debug_QGSM 
    G4QGSMSplitableHadron * at =(G4QGSMSplitableHadron*)theProjectile;
    std::deque<G4Parton *>color=at->GetColorPartons();
    std::deque<G4Parton *>anticolor=at->GetAntiColorPartons();
    G4cout << " proj. color/anti size " << color.size() << " / " << anticolor.size() << G4endl;
    std::deque<G4Parton *>::iterator p_iter;
    G4LorentzVector colmom(0.,0.,0.,0.);
    for ( p_iter=color.begin(); p_iter!= color.end(); ++p_iter){
      G4cout << "proj color : "<<(*p_iter)->GetPDGcode() << ", mom= "<< (*p_iter)->Get4Momentum()<<G4endl;
      colmom+=(*p_iter)->Get4Momentum();
    }
    G4LorentzVector anticolmom(0.,0.,0.,0.);
    for ( p_iter=anticolor.begin(); p_iter!= anticolor.end(); ++p_iter){
      G4cout << "proj antic : "<<(*p_iter)->GetPDGcode() << ", mom= "<< (*p_iter)->Get4Momentum()<<G4endl;
      anticolmom+=(*p_iter)->Get4Momentum();
    }
    G4cout << " proj. color/anti mom " << colmom << " / " << anticolmom  << " Sum: " << colmom+anticolmom <<G4endl;
    #endif
  } else {	 
    G4cout << " none " << G4endl;
  }	    

  if ( theTarget ) {
    G4cout <<  "Target     " << theTarget->GetDefinition()->GetPDGEncoding()
	   << "  " << theTarget->Get4Momentum()<< G4endl;
    mom+=theTarget->Get4Momentum();
    #ifdef debug_QGSM 
    G4QGSMSplitableHadron * at =(G4QGSMSplitableHadron*)theTarget;
    std::deque<G4Parton *>color=at->GetColorPartons();
    std::deque<G4Parton *>anticolor=at->GetAntiColorPartons();
    G4cout << " target color/anti size " << color.size() << " / " << anticolor.size() << G4endl;
    std::deque<G4Parton *>::iterator p_iter;
    G4LorentzVector colmom(0.,0.,0.,0.);
    for ( p_iter=color.begin(); p_iter!= color.end(); ++p_iter){
      G4cout << "target color : "<<(*p_iter)->GetPDGcode() << ", mom= "<< (*p_iter)->Get4Momentum()<<G4endl;
      colmom+=(*p_iter)->Get4Momentum();
    }
    G4LorentzVector anticolmom(0.,0.,0.,0.);
    for ( p_iter=anticolor.begin(); p_iter!= anticolor.end(); ++p_iter){
      G4cout << "target antic : "<<(*p_iter)->GetPDGcode() << ", mom= "<< (*p_iter)->Get4Momentum()<<G4endl;
      anticolmom+=(*p_iter)->Get4Momentum();
    }
    G4cout << " target color/anti mom " << colmom << " / " << anticolmom  << " Sum: " << colmom+anticolmom <<G4endl;
    #endif
  } else {
    G4cout << " none " << G4endl;
  }	    
  G4cout << "total 4-mom of interaction content " << mom << G4endl;
}      
#endif      

#endif

