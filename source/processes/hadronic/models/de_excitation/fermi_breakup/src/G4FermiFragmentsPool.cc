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
//


#include "G4FermiFragmentsPool.hh"

G4bool G4FermiFragmentsPool::MapIsEmpty(true);


std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , std::less<const std::pair<G4int,G4int> > >  &
G4FermiFragmentsPool::GetMap()
{
  //                                                             A  Z  Pol  ExcitE
  static const G4StableFermiFragment Fragment00(  1, 0,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment01(  1, 1,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment02(  2, 1,  3,  0.00*keV );
  static const G4StableFermiFragment Fragment03(  3, 1,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment04(  3, 2,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment05(  4, 2,  1,  0.00*keV );
  static const G4He5FermiFragment    Fragment06(  5, 2,  4, 16.76*keV ); // He5
  static const G4Li5FermiFragment    Fragment07(  5, 3,  4, 16.66*keV ); // Li5
  static const G4StableFermiFragment Fragment08(  6, 2,  1,  0.00*keV );
  static const G4StableFermiFragment Fragment09(  6, 3,  3,  0.00*keV );
  
  static const G4StableFermiFragment Fragment10(  6, 3,  1,  3.56*keV );
  static const G4StableFermiFragment Fragment11(  7, 3,  4,  0.00*keV );
  static const G4StableFermiFragment Fragment12(  7, 3,  2,  0.48*keV );
  static const G4StableFermiFragment Fragment13(  7, 4,  4,  0.00*keV );
  static const G4StableFermiFragment Fragment14(  7, 4,  2,  0.43*keV );
  static const G4StableFermiFragment Fragment15(  8, 3,  5,  0.00*keV );
  static const G4StableFermiFragment Fragment16(  8, 3,  3,  0.98*keV );
  static const G4Be8FermiFragment    Fragment17(  8, 4,  1,  0.00*keV ); // Be8
  static const G4StableFermiFragment Fragment18(  9, 4,  4,  0.00*keV );
  static const G4B9FermiFragment     Fragment19(  9, 5,  4,  0.00*keV ); // B9
  
  static const G4StableFermiFragment Fragment20( 10, 4,  1,  0.00*keV );
  static const G4StableFermiFragment Fragment21( 10, 4,  5,  3.37*keV );
  static const G4StableFermiFragment Fragment22( 10, 4,  8,  5.96*keV );
  static const G4StableFermiFragment Fragment23( 10, 4,  1,  6.18*keV );
  static const G4StableFermiFragment Fragment24( 10, 4,  5,  6.26*keV );
  static const G4StableFermiFragment Fragment25( 10, 5,  7,  0.00*keV );
  static const G4StableFermiFragment Fragment26( 10, 5,  3,  0.72*keV );
  static const G4StableFermiFragment Fragment27( 10, 5,  1,  1.74*keV );
  static const G4StableFermiFragment Fragment28( 10, 5,  3,  2.15*keV );
  static const G4StableFermiFragment Fragment29( 10, 5,  5,  3.59*keV );
  
  static const G4StableFermiFragment Fragment30( 10, 6,  3,  0.00*keV );
  static const G4StableFermiFragment Fragment31( 10, 6,  5,  3.35*keV );
  static const G4StableFermiFragment Fragment32( 11, 5,  4,  0.00*keV );
  static const G4StableFermiFragment Fragment33( 11, 5,  2,  2.13*keV );
  static const G4StableFermiFragment Fragment34( 11, 5,  6,  4.44*keV );
  static const G4StableFermiFragment Fragment35( 11, 5,  4,  5.02*keV );
  static const G4StableFermiFragment Fragment36( 11, 5, 10,  6.76*keV );
  static const G4StableFermiFragment Fragment37( 11, 5,  6,  7.29*keV );
  static const G4StableFermiFragment Fragment38( 11, 5,  4,  7.98*keV );
  static const G4StableFermiFragment Fragment39( 11, 5,  6,  8.56*keV );
  
  static const G4StableFermiFragment Fragment40( 11, 6,  4,  0.00*keV );
  static const G4StableFermiFragment Fragment41( 11, 6,  2,  2.00*keV );
  static const G4StableFermiFragment Fragment42( 11, 6,  6,  4.32*keV );
  static const G4StableFermiFragment Fragment43( 11, 6,  4,  4.80*keV );
  static const G4StableFermiFragment Fragment44( 11, 6,  2,  6.34*keV );
  static const G4StableFermiFragment Fragment45( 11, 6,  8,  6.48*keV );
  static const G4StableFermiFragment Fragment46( 11, 6,  6,  6.90*keV );
  static const G4StableFermiFragment Fragment47( 11, 6,  4,  7.50*keV );
  static const G4StableFermiFragment Fragment48( 11, 6,  4,  8.10*keV );
  static const G4StableFermiFragment Fragment49( 11, 6,  6,  8.42*keV );
  
  static const G4StableFermiFragment Fragment50( 11, 6,  8,  8.66*keV );
  static const G4StableFermiFragment Fragment51( 12, 5,  3,  0.00*keV );
  static const G4StableFermiFragment Fragment52( 12, 5,  5,  0.95*keV );
  static const G4StableFermiFragment Fragment53( 12, 5,  5,  1.67*keV );
  static const G4StableFermiFragment Fragment54( 12, 5,  4,  2.65*keV );
  static const G4StableFermiFragment Fragment55( 12, 6,  1,  0.00*keV );
  static const G4StableFermiFragment Fragment56( 12, 6,  5,  4.44*keV );
  static const G4StableFermiFragment Fragment57( 13, 6,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment58( 13, 6,  2,  3.09*keV );
  static const G4StableFermiFragment Fragment59( 13, 6,  4,  3.68*keV );
  
  static const G4StableFermiFragment Fragment60( 13, 6,  6,  3.85*keV );
  static const G4StableFermiFragment Fragment61( 13, 7,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment62( 14, 6,  1,  0.00*keV );
  static const G4StableFermiFragment Fragment63( 14, 6,  3,  6.09*keV );
  static const G4StableFermiFragment Fragment64( 14, 6,  8,  6.69*keV );
  static const G4StableFermiFragment Fragment65( 14, 6,  6,  6.96*keV );
  static const G4StableFermiFragment Fragment66( 14, 6,  5,  7.34*keV );
  static const G4StableFermiFragment Fragment67( 14, 7,  3,  0.00*keV );
  static const G4StableFermiFragment Fragment68( 14, 7,  1,  2.31*keV );
  static const G4StableFermiFragment Fragment69( 14, 7,  3,  3.95*keV );
  
  static const G4StableFermiFragment Fragment70( 14, 7,  1,  4.92*keV );
  static const G4StableFermiFragment Fragment71( 14, 7,  5,  5.11*keV );
  static const G4StableFermiFragment Fragment72( 14, 7,  3,  5.69*keV );
  static const G4StableFermiFragment Fragment73( 14, 7,  7,  5.83*keV );
  static const G4StableFermiFragment Fragment74( 14, 7,  3,  6.20*keV );
  static const G4StableFermiFragment Fragment75( 14, 7,  7,  6.44*keV );
  static const G4StableFermiFragment Fragment76( 14, 7,  5,  7.03*keV );
  static const G4StableFermiFragment Fragment77( 15, 7,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment78( 15, 7,  8,  5.28*keV );
  static const G4StableFermiFragment Fragment79( 15, 7,  4,  6.32*keV );
  
  static const G4StableFermiFragment Fragment80( 15, 7, 10,  7.22*keV );
  static const G4StableFermiFragment Fragment81( 15, 7,  8,  7.57*keV );
  static const G4StableFermiFragment Fragment82( 15, 7,  2,  8.31*keV );
  static const G4StableFermiFragment Fragment83( 15, 7,  4,  8.57*keV );
  static const G4StableFermiFragment Fragment84( 15, 7, 14,  9.15*keV );
  static const G4StableFermiFragment Fragment85( 15, 7, 14,  9.79*keV );
  static const G4StableFermiFragment Fragment86( 15, 7,  8, 10.00*keV );
  static const G4StableFermiFragment Fragment87( 15, 8,  2,  0.00*keV );
  static const G4StableFermiFragment Fragment88( 15, 8,  8,  5.22*keV );
  static const G4StableFermiFragment Fragment89( 15, 8,  4,  6.18*keV );
  
  static const G4StableFermiFragment Fragment90( 15, 8, 10,  6.83*keV );
  static const G4StableFermiFragment Fragment91( 15, 8,  8,  7.28*keV );
  static const G4StableFermiFragment Fragment92( 16, 7,  5,  0.00*keV );
  static const G4StableFermiFragment Fragment93( 16, 7,  1,  0.12*keV );
  static const G4StableFermiFragment Fragment94( 16, 7,  7,  0.30*keV );
  static const G4StableFermiFragment Fragment95( 16, 7,  3,  0.40*keV );
  static const G4StableFermiFragment Fragment96( 16, 8,  1,  0.00*keV );
  static const G4StableFermiFragment Fragment97( 16, 8,  8,  6.10*keV );
  static const G4StableFermiFragment Fragment98( 16, 8,  5,  6.92*keV );
  static const G4StableFermiFragment Fragment99( 16, 8,  3,  7.12*keV );

  static std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , std::less<const std::pair<G4int,G4int> > >  
    theMapOfFragments;

  if (MapIsEmpty) 
    {
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment00.GetA(),Fragment00.GetZ()),&Fragment00));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment01.GetA(),Fragment01.GetZ()),&Fragment01));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment02.GetA(),Fragment02.GetZ()),&Fragment02));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment03.GetA(),Fragment03.GetZ()),&Fragment03));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment04.GetA(),Fragment04.GetZ()),&Fragment04));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment05.GetA(),Fragment05.GetZ()),&Fragment05));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment06.GetA(),Fragment06.GetZ()),&Fragment06));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment07.GetA(),Fragment07.GetZ()),&Fragment07));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment08.GetA(),Fragment08.GetZ()),&Fragment08));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment09.GetA(),Fragment09.GetZ()),&Fragment09));

      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment10.GetA(),Fragment10.GetZ()),&Fragment10));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment11.GetA(),Fragment11.GetZ()),&Fragment11));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment12.GetA(),Fragment12.GetZ()),&Fragment12));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment13.GetA(),Fragment13.GetZ()),&Fragment13));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment14.GetA(),Fragment14.GetZ()),&Fragment14));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment15.GetA(),Fragment15.GetZ()),&Fragment15));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment16.GetA(),Fragment16.GetZ()),&Fragment16));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment17.GetA(),Fragment17.GetZ()),&Fragment17));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment18.GetA(),Fragment18.GetZ()),&Fragment18));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment19.GetA(),Fragment19.GetZ()),&Fragment19));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment20.GetA(),Fragment20.GetZ()),&Fragment20));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment21.GetA(),Fragment21.GetZ()),&Fragment21));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment22.GetA(),Fragment22.GetZ()),&Fragment22));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment23.GetA(),Fragment23.GetZ()),&Fragment23));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment24.GetA(),Fragment24.GetZ()),&Fragment24));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment25.GetA(),Fragment25.GetZ()),&Fragment25));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment26.GetA(),Fragment26.GetZ()),&Fragment26));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment27.GetA(),Fragment27.GetZ()),&Fragment27));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment28.GetA(),Fragment28.GetZ()),&Fragment28));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment29.GetA(),Fragment29.GetZ()),&Fragment29));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment30.GetA(),Fragment30.GetZ()),&Fragment30));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment31.GetA(),Fragment31.GetZ()),&Fragment31));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment32.GetA(),Fragment32.GetZ()),&Fragment32));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment33.GetA(),Fragment33.GetZ()),&Fragment33));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment34.GetA(),Fragment34.GetZ()),&Fragment34));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment35.GetA(),Fragment35.GetZ()),&Fragment35));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment36.GetA(),Fragment36.GetZ()),&Fragment36));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment37.GetA(),Fragment37.GetZ()),&Fragment37));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment38.GetA(),Fragment38.GetZ()),&Fragment38));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment39.GetA(),Fragment39.GetZ()),&Fragment39));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment40.GetA(),Fragment40.GetZ()),&Fragment40));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment41.GetA(),Fragment41.GetZ()),&Fragment41));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment42.GetA(),Fragment42.GetZ()),&Fragment42));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment43.GetA(),Fragment43.GetZ()),&Fragment43));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment44.GetA(),Fragment44.GetZ()),&Fragment44));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment45.GetA(),Fragment45.GetZ()),&Fragment45));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment46.GetA(),Fragment46.GetZ()),&Fragment46));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment47.GetA(),Fragment47.GetZ()),&Fragment47));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment48.GetA(),Fragment48.GetZ()),&Fragment48));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment49.GetA(),Fragment49.GetZ()),&Fragment49));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment50.GetA(),Fragment50.GetZ()),&Fragment50));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment51.GetA(),Fragment51.GetZ()),&Fragment51));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment52.GetA(),Fragment52.GetZ()),&Fragment52));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment53.GetA(),Fragment53.GetZ()),&Fragment53));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment54.GetA(),Fragment54.GetZ()),&Fragment54));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment55.GetA(),Fragment55.GetZ()),&Fragment55));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment56.GetA(),Fragment56.GetZ()),&Fragment56));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment57.GetA(),Fragment57.GetZ()),&Fragment57));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment58.GetA(),Fragment58.GetZ()),&Fragment58));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment59.GetA(),Fragment59.GetZ()),&Fragment59));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment60.GetA(),Fragment60.GetZ()),&Fragment60));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment61.GetA(),Fragment61.GetZ()),&Fragment61));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment62.GetA(),Fragment62.GetZ()),&Fragment62));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment63.GetA(),Fragment63.GetZ()),&Fragment63));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment64.GetA(),Fragment64.GetZ()),&Fragment64));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment65.GetA(),Fragment65.GetZ()),&Fragment65));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment66.GetA(),Fragment66.GetZ()),&Fragment66));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment67.GetA(),Fragment67.GetZ()),&Fragment67));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment68.GetA(),Fragment68.GetZ()),&Fragment68));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment69.GetA(),Fragment69.GetZ()),&Fragment69));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment70.GetA(),Fragment70.GetZ()),&Fragment70));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment71.GetA(),Fragment71.GetZ()),&Fragment71));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment72.GetA(),Fragment72.GetZ()),&Fragment72));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment73.GetA(),Fragment73.GetZ()),&Fragment73));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment74.GetA(),Fragment74.GetZ()),&Fragment74));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment75.GetA(),Fragment75.GetZ()),&Fragment75));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment76.GetA(),Fragment76.GetZ()),&Fragment76));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment77.GetA(),Fragment77.GetZ()),&Fragment77));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment78.GetA(),Fragment78.GetZ()),&Fragment78));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment79.GetA(),Fragment79.GetZ()),&Fragment79));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment80.GetA(),Fragment80.GetZ()),&Fragment80));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment81.GetA(),Fragment81.GetZ()),&Fragment81));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment82.GetA(),Fragment82.GetZ()),&Fragment82));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment83.GetA(),Fragment83.GetZ()),&Fragment83));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment84.GetA(),Fragment84.GetZ()),&Fragment84));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment85.GetA(),Fragment85.GetZ()),&Fragment85));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment86.GetA(),Fragment86.GetZ()),&Fragment86));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment87.GetA(),Fragment87.GetZ()),&Fragment87));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment88.GetA(),Fragment88.GetZ()),&Fragment88));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment89.GetA(),Fragment89.GetZ()),&Fragment89));
      
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment90.GetA(),Fragment90.GetZ()),&Fragment90));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment91.GetA(),Fragment91.GetZ()),&Fragment91));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment92.GetA(),Fragment92.GetZ()),&Fragment92));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment93.GetA(),Fragment93.GetZ()),&Fragment93));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment94.GetA(),Fragment94.GetZ()),&Fragment94));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment95.GetA(),Fragment95.GetZ()),&Fragment95));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment96.GetA(),Fragment96.GetZ()),&Fragment96));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment97.GetA(),Fragment97.GetZ()),&Fragment97));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment98.GetA(),Fragment98.GetZ()),&Fragment98));
      theMapOfFragments.insert(std::pair<const std::pair<G4int,G4int>, const G4VFermiFragment* >(std::pair<G4int,G4int>(Fragment99.GetA(),Fragment99.GetZ()),&Fragment99));
      MapIsEmpty = false;
    }

  return theMapOfFragments;
}


G4FermiFragmentsPool::G4FermiFragmentsPool()
{
}



G4FermiFragmentsPool::~G4FermiFragmentsPool()
{
}





G4FermiFragmentsPool::G4FermiFragmentsPool(const G4FermiFragmentsPool&)
{
  // It is meant to not be accesable
}

const G4FermiFragmentsPool & G4FermiFragmentsPool::operator=(const G4FermiFragmentsPool& )
{
  // It is meant to not be accesable
  return *this;
}

G4bool G4FermiFragmentsPool::operator==(const G4FermiFragmentsPool&) const
{
  // It is meant to not be accesable
  return false;
}

G4bool G4FermiFragmentsPool::operator!=(const G4FermiFragmentsPool&) const
{
  // It is meant to not be accesable
  return true;
}
