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
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4FermiFragmentsPool_hh 
#define G4FermiFragmentsPool_hh

#include "G4VFermiFragment.hh"
#include "G4StableFermiFragment.hh"
#include "G4B9FermiFragment.hh"
#include "G4Be8FermiFragment.hh"
#include "G4He5FermiFragment.hh"
#include "G4Li5FermiFragment.hh"

#include <map>
#include <functional>

class G4FermiFragmentsPool
{
public:
  G4FermiFragmentsPool();
  ~G4FermiFragmentsPool();


  static G4int Count(const std::pair<G4int,G4int>& az) 
  {
    return theMapOfFragments.count(az);
  }
  
  static std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
		       std::less<std::pair<G4int,G4int> > >::iterator 
  LowerBound(const std::pair<G4int,G4int> & az) 
  {
    return theMapOfFragments.lower_bound(az);
  }

  static std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
			 std::less<std::pair<G4int,G4int> > >::iterator 
  UpperBound(const std::pair<G4int,G4int> & az) 
  {
    return theMapOfFragments.upper_bound(az);
  }

  


  
private:
  G4FermiFragmentsPool(const G4FermiFragmentsPool&);
  const G4FermiFragmentsPool & operator=(const G4FermiFragmentsPool&);
  G4bool operator==(const G4FermiFragmentsPool&) const;
  G4bool operator!=(const G4FermiFragmentsPool&) const;

private:
  static const G4StableFermiFragment Fragment00;
  static const G4StableFermiFragment Fragment01;
  static const G4StableFermiFragment Fragment02;
  static const G4StableFermiFragment Fragment03;
  static const G4StableFermiFragment Fragment04;
  static const G4StableFermiFragment Fragment05;
  static const G4He5FermiFragment    Fragment06;   // He5
  static const G4Li5FermiFragment    Fragment07;   // Li5
  static const G4StableFermiFragment Fragment08; 
  static const G4StableFermiFragment Fragment09;

  static const G4StableFermiFragment Fragment10;
  static const G4StableFermiFragment Fragment11;
  static const G4StableFermiFragment Fragment12;
  static const G4StableFermiFragment Fragment13;
  static const G4StableFermiFragment Fragment14;
  static const G4StableFermiFragment Fragment15;
  static const G4StableFermiFragment Fragment16;
  static const G4Be8FermiFragment    Fragment17;  // Be8
  static const G4StableFermiFragment Fragment18;
  static const G4B9FermiFragment     Fragment19;  // B9

  static const G4StableFermiFragment Fragment20;
  static const G4StableFermiFragment Fragment21;
  static const G4StableFermiFragment Fragment22;
  static const G4StableFermiFragment Fragment23;
  static const G4StableFermiFragment Fragment24;
  static const G4StableFermiFragment Fragment25;
  static const G4StableFermiFragment Fragment26;
  static const G4StableFermiFragment Fragment27;
  static const G4StableFermiFragment Fragment28;
  static const G4StableFermiFragment Fragment29;

  static const G4StableFermiFragment Fragment30;
  static const G4StableFermiFragment Fragment31;
  static const G4StableFermiFragment Fragment32;
  static const G4StableFermiFragment Fragment33;
  static const G4StableFermiFragment Fragment34;
  static const G4StableFermiFragment Fragment35;
  static const G4StableFermiFragment Fragment36;
  static const G4StableFermiFragment Fragment37;
  static const G4StableFermiFragment Fragment38;
  static const G4StableFermiFragment Fragment39;

  static const G4StableFermiFragment Fragment40;
  static const G4StableFermiFragment Fragment41;
  static const G4StableFermiFragment Fragment42;
  static const G4StableFermiFragment Fragment43;
  static const G4StableFermiFragment Fragment44;
  static const G4StableFermiFragment Fragment45;
  static const G4StableFermiFragment Fragment46;
  static const G4StableFermiFragment Fragment47;
  static const G4StableFermiFragment Fragment48;
  static const G4StableFermiFragment Fragment49;

  static const G4StableFermiFragment Fragment50;
  static const G4StableFermiFragment Fragment51;
  static const G4StableFermiFragment Fragment52;
  static const G4StableFermiFragment Fragment53;
  static const G4StableFermiFragment Fragment54;
  static const G4StableFermiFragment Fragment55;
  static const G4StableFermiFragment Fragment56;
  static const G4StableFermiFragment Fragment57;
  static const G4StableFermiFragment Fragment58;
  static const G4StableFermiFragment Fragment59;

  static const G4StableFermiFragment Fragment60;
  static const G4StableFermiFragment Fragment61;
  static const G4StableFermiFragment Fragment62;
  static const G4StableFermiFragment Fragment63;
  static const G4StableFermiFragment Fragment64;
  static const G4StableFermiFragment Fragment65;
  static const G4StableFermiFragment Fragment66;
  static const G4StableFermiFragment Fragment67;
  static const G4StableFermiFragment Fragment68;
  static const G4StableFermiFragment Fragment69;

  static const G4StableFermiFragment Fragment70;
  static const G4StableFermiFragment Fragment71;
  static const G4StableFermiFragment Fragment72;
  static const G4StableFermiFragment Fragment73;
  static const G4StableFermiFragment Fragment74;
  static const G4StableFermiFragment Fragment75;
  static const G4StableFermiFragment Fragment76;
  static const G4StableFermiFragment Fragment77;
  static const G4StableFermiFragment Fragment78;
  static const G4StableFermiFragment Fragment79;

  static const G4StableFermiFragment Fragment80;
  static const G4StableFermiFragment Fragment81;
  static const G4StableFermiFragment Fragment82;
  static const G4StableFermiFragment Fragment83;
  static const G4StableFermiFragment Fragment84;
  static const G4StableFermiFragment Fragment85;
  static const G4StableFermiFragment Fragment86;
  static const G4StableFermiFragment Fragment87;
  static const G4StableFermiFragment Fragment88;
  static const G4StableFermiFragment Fragment89;

  static const G4StableFermiFragment Fragment90;
  static const G4StableFermiFragment Fragment91;
  static const G4StableFermiFragment Fragment92;
  static const G4StableFermiFragment Fragment93;
  static const G4StableFermiFragment Fragment94;
  static const G4StableFermiFragment Fragment95;
  static const G4StableFermiFragment Fragment96;
  static const G4StableFermiFragment Fragment97;
  static const G4StableFermiFragment Fragment98;
  static const G4StableFermiFragment Fragment99;


  static std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , 
		       std::less<const std::pair<G4int,G4int> > >  theMapOfFragments;
};
#endif

