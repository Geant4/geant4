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

//                                                             A  Z  Pol  ExcitE
const G4StableFermiFragment G4FermiFragmentsPool::Fragment00(  1, 0,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment01(  1, 1,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment02(  2, 1,  3,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment03(  3, 1,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment04(  3, 2,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment05(  4, 2,  1,  0.00*keV );
const G4He5FermiFragment    G4FermiFragmentsPool::Fragment06(  5, 2,  4, 16.76*keV ); // He5
const G4Li5FermiFragment    G4FermiFragmentsPool::Fragment07(  5, 3,  4, 16.66*keV ); // Li5
const G4StableFermiFragment G4FermiFragmentsPool::Fragment08(  6, 2,  1,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment09(  6, 3,  3,  0.00*keV );
			                    
const G4StableFermiFragment G4FermiFragmentsPool::Fragment10(  6, 3,  1,  3.56*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment11(  7, 3,  4,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment12(  7, 3,  2,  0.48*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment13(  7, 4,  4,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment14(  7, 4,  2,  0.43*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment15(  8, 3,  5,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment16(  8, 3,  3,  0.98*keV );
const G4Be8FermiFragment    G4FermiFragmentsPool::Fragment17(  8, 4,  1,  0.00*keV ); // Be8
const G4StableFermiFragment G4FermiFragmentsPool::Fragment18(  9, 4,  4,  0.00*keV );
const G4B9FermiFragment     G4FermiFragmentsPool::Fragment19(  9, 5,  4,  0.00*keV ); // B9
			                     
const G4StableFermiFragment G4FermiFragmentsPool::Fragment20( 10, 4,  1,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment21( 10, 4,  5,  3.37*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment22( 10, 4,  8,  5.96*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment23( 10, 4,  1,  6.18*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment24( 10, 4,  5,  6.26*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment25( 10, 5,  7,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment26( 10, 5,  3,  0.72*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment27( 10, 5,  1,  1.74*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment28( 10, 5,  3,  2.15*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment29( 10, 5,  5,  3.59*keV );
			                     
const G4StableFermiFragment G4FermiFragmentsPool::Fragment30( 10, 6,  3,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment31( 10, 6,  5,  3.35*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment32( 11, 5,  4,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment33( 11, 5,  2,  2.13*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment34( 11, 5,  6,  4.44*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment35( 11, 5,  4,  5.02*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment36( 11, 5, 10,  6.76*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment37( 11, 5,  6,  7.29*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment38( 11, 5,  4,  7.98*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment39( 11, 5,  6,  8.56*keV );
        					   
const G4StableFermiFragment G4FermiFragmentsPool::Fragment40( 11, 6,  4,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment41( 11, 6,  2,  2.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment42( 11, 6,  6,  4.32*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment43( 11, 6,  4,  4.80*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment44( 11, 6,  2,  6.34*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment45( 11, 6,  8,  6.48*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment46( 11, 6,  6,  6.90*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment47( 11, 6,  4,  7.50*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment48( 11, 6,  4,  8.10*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment49( 11, 6,  6,  8.42*keV );
        					   
const G4StableFermiFragment G4FermiFragmentsPool::Fragment50( 11, 6,  8,  8.66*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment51( 12, 5,  3,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment52( 12, 5,  5,  0.95*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment53( 12, 5,  5,  1.67*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment54( 12, 5,  4,  2.65*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment55( 12, 6,  1,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment56( 12, 6,  5,  4.44*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment57( 13, 6,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment58( 13, 6,  2,  3.09*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment59( 13, 6,  4,  3.68*keV );
			                     
const G4StableFermiFragment G4FermiFragmentsPool::Fragment60( 13, 6,  6,  3.85*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment61( 13, 7,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment62( 14, 6,  1,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment63( 14, 6,  3,  6.09*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment64( 14, 6,  8,  6.69*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment65( 14, 6,  6,  6.96*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment66( 14, 6,  5,  7.34*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment67( 14, 7,  3,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment68( 14, 7,  1,  2.31*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment69( 14, 7,  3,  3.95*keV );
        					   
const G4StableFermiFragment G4FermiFragmentsPool::Fragment70( 14, 7,  1,  4.92*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment71( 14, 7,  5,  5.11*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment72( 14, 7,  3,  5.69*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment73( 14, 7,  7,  5.83*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment74( 14, 7,  3,  6.20*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment75( 14, 7,  7,  6.44*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment76( 14, 7,  5,  7.03*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment77( 15, 7,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment78( 15, 7,  8,  5.28*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment79( 15, 7,  4,  6.32*keV );
        					   
const G4StableFermiFragment G4FermiFragmentsPool::Fragment80( 15, 7, 10,  7.22*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment81( 15, 7,  8,  7.57*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment82( 15, 7,  2,  8.31*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment83( 15, 7,  4,  8.57*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment84( 15, 7, 14,  9.15*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment85( 15, 7, 14,  9.79*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment86( 15, 7,  8, 10.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment87( 15, 8,  2,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment88( 15, 8,  8,  5.22*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment89( 15, 8,  4,  6.18*keV );
			                     
const G4StableFermiFragment G4FermiFragmentsPool::Fragment90( 15, 8, 10,  6.83*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment91( 15, 8,  8,  7.28*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment92( 16, 7,  5,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment93( 16, 7,  1,  0.12*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment94( 16, 7,  7,  0.30*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment95( 16, 7,  3,  0.40*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment96( 16, 8,  1,  0.00*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment97( 16, 8,  8,  6.10*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment98( 16, 8,  5,  6.92*keV );
const G4StableFermiFragment G4FermiFragmentsPool::Fragment99( 16, 8,  3,  7.12*keV );


std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , 
	      std::less<const std::pair<G4int,G4int> > >  G4FermiFragmentsPool::theMapOfFragments;


G4FermiFragmentsPool::G4FermiFragmentsPool()
{
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment00.GetA(),G4FermiFragmentsPool::Fragment00.GetZ()),&G4FermiFragmentsPool::Fragment00));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment01.GetA(),G4FermiFragmentsPool::Fragment01.GetZ()),&G4FermiFragmentsPool::Fragment01));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment02.GetA(),G4FermiFragmentsPool::Fragment02.GetZ()),&G4FermiFragmentsPool::Fragment02));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment03.GetA(),G4FermiFragmentsPool::Fragment03.GetZ()),&G4FermiFragmentsPool::Fragment03));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment04.GetA(),G4FermiFragmentsPool::Fragment04.GetZ()),&G4FermiFragmentsPool::Fragment04));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment05.GetA(),G4FermiFragmentsPool::Fragment05.GetZ()),&G4FermiFragmentsPool::Fragment05));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment06.GetA(),G4FermiFragmentsPool::Fragment06.GetZ()),&G4FermiFragmentsPool::Fragment06));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment07.GetA(),G4FermiFragmentsPool::Fragment07.GetZ()),&G4FermiFragmentsPool::Fragment07));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment08.GetA(),G4FermiFragmentsPool::Fragment08.GetZ()),&G4FermiFragmentsPool::Fragment08));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment09.GetA(),G4FermiFragmentsPool::Fragment09.GetZ()),&G4FermiFragmentsPool::Fragment09));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment10.GetA(),G4FermiFragmentsPool::Fragment10.GetZ()),&G4FermiFragmentsPool::Fragment10));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment11.GetA(),G4FermiFragmentsPool::Fragment11.GetZ()),&G4FermiFragmentsPool::Fragment11));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment12.GetA(),G4FermiFragmentsPool::Fragment12.GetZ()),&G4FermiFragmentsPool::Fragment12));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment13.GetA(),G4FermiFragmentsPool::Fragment13.GetZ()),&G4FermiFragmentsPool::Fragment13));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment14.GetA(),G4FermiFragmentsPool::Fragment14.GetZ()),&G4FermiFragmentsPool::Fragment14));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment15.GetA(),G4FermiFragmentsPool::Fragment15.GetZ()),&G4FermiFragmentsPool::Fragment15));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment16.GetA(),G4FermiFragmentsPool::Fragment16.GetZ()),&G4FermiFragmentsPool::Fragment16));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment17.GetA(),G4FermiFragmentsPool::Fragment17.GetZ()),&G4FermiFragmentsPool::Fragment17));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment18.GetA(),G4FermiFragmentsPool::Fragment18.GetZ()),&G4FermiFragmentsPool::Fragment18));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment19.GetA(),G4FermiFragmentsPool::Fragment19.GetZ()),&G4FermiFragmentsPool::Fragment19));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment20.GetA(),G4FermiFragmentsPool::Fragment20.GetZ()),&G4FermiFragmentsPool::Fragment20));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment21.GetA(),G4FermiFragmentsPool::Fragment21.GetZ()),&G4FermiFragmentsPool::Fragment21));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment22.GetA(),G4FermiFragmentsPool::Fragment22.GetZ()),&G4FermiFragmentsPool::Fragment22));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment23.GetA(),G4FermiFragmentsPool::Fragment23.GetZ()),&G4FermiFragmentsPool::Fragment23));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment24.GetA(),G4FermiFragmentsPool::Fragment24.GetZ()),&G4FermiFragmentsPool::Fragment24));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment25.GetA(),G4FermiFragmentsPool::Fragment25.GetZ()),&G4FermiFragmentsPool::Fragment25));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment26.GetA(),G4FermiFragmentsPool::Fragment26.GetZ()),&G4FermiFragmentsPool::Fragment26));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment27.GetA(),G4FermiFragmentsPool::Fragment27.GetZ()),&G4FermiFragmentsPool::Fragment27));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment28.GetA(),G4FermiFragmentsPool::Fragment28.GetZ()),&G4FermiFragmentsPool::Fragment28));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment29.GetA(),G4FermiFragmentsPool::Fragment29.GetZ()),&G4FermiFragmentsPool::Fragment29));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment30.GetA(),G4FermiFragmentsPool::Fragment30.GetZ()),&G4FermiFragmentsPool::Fragment30));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment31.GetA(),G4FermiFragmentsPool::Fragment31.GetZ()),&G4FermiFragmentsPool::Fragment31));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment32.GetA(),G4FermiFragmentsPool::Fragment32.GetZ()),&G4FermiFragmentsPool::Fragment32));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment33.GetA(),G4FermiFragmentsPool::Fragment33.GetZ()),&G4FermiFragmentsPool::Fragment33));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment34.GetA(),G4FermiFragmentsPool::Fragment34.GetZ()),&G4FermiFragmentsPool::Fragment34));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment35.GetA(),G4FermiFragmentsPool::Fragment35.GetZ()),&G4FermiFragmentsPool::Fragment35));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment36.GetA(),G4FermiFragmentsPool::Fragment36.GetZ()),&G4FermiFragmentsPool::Fragment36));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment37.GetA(),G4FermiFragmentsPool::Fragment37.GetZ()),&G4FermiFragmentsPool::Fragment37));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment38.GetA(),G4FermiFragmentsPool::Fragment38.GetZ()),&G4FermiFragmentsPool::Fragment38));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment39.GetA(),G4FermiFragmentsPool::Fragment39.GetZ()),&G4FermiFragmentsPool::Fragment39));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment40.GetA(),G4FermiFragmentsPool::Fragment40.GetZ()),&G4FermiFragmentsPool::Fragment40));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment41.GetA(),G4FermiFragmentsPool::Fragment41.GetZ()),&G4FermiFragmentsPool::Fragment41));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment42.GetA(),G4FermiFragmentsPool::Fragment42.GetZ()),&G4FermiFragmentsPool::Fragment42));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment43.GetA(),G4FermiFragmentsPool::Fragment43.GetZ()),&G4FermiFragmentsPool::Fragment43));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment44.GetA(),G4FermiFragmentsPool::Fragment44.GetZ()),&G4FermiFragmentsPool::Fragment44));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment45.GetA(),G4FermiFragmentsPool::Fragment45.GetZ()),&G4FermiFragmentsPool::Fragment45));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment46.GetA(),G4FermiFragmentsPool::Fragment46.GetZ()),&G4FermiFragmentsPool::Fragment46));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment47.GetA(),G4FermiFragmentsPool::Fragment47.GetZ()),&G4FermiFragmentsPool::Fragment47));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment48.GetA(),G4FermiFragmentsPool::Fragment48.GetZ()),&G4FermiFragmentsPool::Fragment48));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment49.GetA(),G4FermiFragmentsPool::Fragment49.GetZ()),&G4FermiFragmentsPool::Fragment49));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment50.GetA(),G4FermiFragmentsPool::Fragment50.GetZ()),&G4FermiFragmentsPool::Fragment50));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment51.GetA(),G4FermiFragmentsPool::Fragment51.GetZ()),&G4FermiFragmentsPool::Fragment51));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment52.GetA(),G4FermiFragmentsPool::Fragment52.GetZ()),&G4FermiFragmentsPool::Fragment52));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment53.GetA(),G4FermiFragmentsPool::Fragment53.GetZ()),&G4FermiFragmentsPool::Fragment53));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment54.GetA(),G4FermiFragmentsPool::Fragment54.GetZ()),&G4FermiFragmentsPool::Fragment54));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment55.GetA(),G4FermiFragmentsPool::Fragment55.GetZ()),&G4FermiFragmentsPool::Fragment55));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment56.GetA(),G4FermiFragmentsPool::Fragment56.GetZ()),&G4FermiFragmentsPool::Fragment56));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment57.GetA(),G4FermiFragmentsPool::Fragment57.GetZ()),&G4FermiFragmentsPool::Fragment57));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment58.GetA(),G4FermiFragmentsPool::Fragment58.GetZ()),&G4FermiFragmentsPool::Fragment58));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment59.GetA(),G4FermiFragmentsPool::Fragment59.GetZ()),&G4FermiFragmentsPool::Fragment59));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment60.GetA(),G4FermiFragmentsPool::Fragment60.GetZ()),&G4FermiFragmentsPool::Fragment60));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment61.GetA(),G4FermiFragmentsPool::Fragment61.GetZ()),&G4FermiFragmentsPool::Fragment61));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment62.GetA(),G4FermiFragmentsPool::Fragment62.GetZ()),&G4FermiFragmentsPool::Fragment62));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment63.GetA(),G4FermiFragmentsPool::Fragment63.GetZ()),&G4FermiFragmentsPool::Fragment63));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment64.GetA(),G4FermiFragmentsPool::Fragment64.GetZ()),&G4FermiFragmentsPool::Fragment64));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment65.GetA(),G4FermiFragmentsPool::Fragment65.GetZ()),&G4FermiFragmentsPool::Fragment65));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment66.GetA(),G4FermiFragmentsPool::Fragment66.GetZ()),&G4FermiFragmentsPool::Fragment66));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment67.GetA(),G4FermiFragmentsPool::Fragment67.GetZ()),&G4FermiFragmentsPool::Fragment67));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment68.GetA(),G4FermiFragmentsPool::Fragment68.GetZ()),&G4FermiFragmentsPool::Fragment68));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment69.GetA(),G4FermiFragmentsPool::Fragment69.GetZ()),&G4FermiFragmentsPool::Fragment69));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment70.GetA(),G4FermiFragmentsPool::Fragment70.GetZ()),&G4FermiFragmentsPool::Fragment70));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment71.GetA(),G4FermiFragmentsPool::Fragment71.GetZ()),&G4FermiFragmentsPool::Fragment71));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment72.GetA(),G4FermiFragmentsPool::Fragment72.GetZ()),&G4FermiFragmentsPool::Fragment72));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment73.GetA(),G4FermiFragmentsPool::Fragment73.GetZ()),&G4FermiFragmentsPool::Fragment73));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment74.GetA(),G4FermiFragmentsPool::Fragment74.GetZ()),&G4FermiFragmentsPool::Fragment74));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment75.GetA(),G4FermiFragmentsPool::Fragment75.GetZ()),&G4FermiFragmentsPool::Fragment75));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment76.GetA(),G4FermiFragmentsPool::Fragment76.GetZ()),&G4FermiFragmentsPool::Fragment76));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment77.GetA(),G4FermiFragmentsPool::Fragment77.GetZ()),&G4FermiFragmentsPool::Fragment77));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment78.GetA(),G4FermiFragmentsPool::Fragment78.GetZ()),&G4FermiFragmentsPool::Fragment78));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment79.GetA(),G4FermiFragmentsPool::Fragment79.GetZ()),&G4FermiFragmentsPool::Fragment79));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment80.GetA(),G4FermiFragmentsPool::Fragment80.GetZ()),&G4FermiFragmentsPool::Fragment80));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment81.GetA(),G4FermiFragmentsPool::Fragment81.GetZ()),&G4FermiFragmentsPool::Fragment81));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment82.GetA(),G4FermiFragmentsPool::Fragment82.GetZ()),&G4FermiFragmentsPool::Fragment82));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment83.GetA(),G4FermiFragmentsPool::Fragment83.GetZ()),&G4FermiFragmentsPool::Fragment83));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment84.GetA(),G4FermiFragmentsPool::Fragment84.GetZ()),&G4FermiFragmentsPool::Fragment84));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment85.GetA(),G4FermiFragmentsPool::Fragment85.GetZ()),&G4FermiFragmentsPool::Fragment85));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment86.GetA(),G4FermiFragmentsPool::Fragment86.GetZ()),&G4FermiFragmentsPool::Fragment86));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment87.GetA(),G4FermiFragmentsPool::Fragment87.GetZ()),&G4FermiFragmentsPool::Fragment87));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment88.GetA(),G4FermiFragmentsPool::Fragment88.GetZ()),&G4FermiFragmentsPool::Fragment88));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment89.GetA(),G4FermiFragmentsPool::Fragment89.GetZ()),&G4FermiFragmentsPool::Fragment89));

  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment90.GetA(),G4FermiFragmentsPool::Fragment90.GetZ()),&G4FermiFragmentsPool::Fragment90));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment91.GetA(),G4FermiFragmentsPool::Fragment91.GetZ()),&G4FermiFragmentsPool::Fragment91));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment92.GetA(),G4FermiFragmentsPool::Fragment92.GetZ()),&G4FermiFragmentsPool::Fragment92));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment93.GetA(),G4FermiFragmentsPool::Fragment93.GetZ()),&G4FermiFragmentsPool::Fragment93));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment94.GetA(),G4FermiFragmentsPool::Fragment94.GetZ()),&G4FermiFragmentsPool::Fragment94));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment95.GetA(),G4FermiFragmentsPool::Fragment95.GetZ()),&G4FermiFragmentsPool::Fragment95));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment96.GetA(),G4FermiFragmentsPool::Fragment96.GetZ()),&G4FermiFragmentsPool::Fragment96));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment97.GetA(),G4FermiFragmentsPool::Fragment97.GetZ()),&G4FermiFragmentsPool::Fragment97));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment98.GetA(),G4FermiFragmentsPool::Fragment98.GetZ()),&G4FermiFragmentsPool::Fragment98));
  theMapOfFragments.insert(std::make_pair(std::make_pair(G4FermiFragmentsPool::Fragment99.GetA(),G4FermiFragmentsPool::Fragment99.GetZ()),&G4FermiFragmentsPool::Fragment99));

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
