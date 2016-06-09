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
// $Id$
//
//
//      ---------------- G4QIsotope header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of the G4QIsotope class of the CHIPS Simulation Branch in GEANT4
// ----------------------------------------------------------------------------
//  Short descriptionIt contains information about natural abundances of stable
//  and long living isotopes and a NEW "Element" can be initialised for any
//  isotope set. Randomization of isotopes of the Natural Elements is hardwired and
//  fast randomization of isotopes of the user defined Elements is a bit slower
//  CrossSectionWeighted randomisation of isotopes is slow (same for Nat and New)
// -------------------------------------------------------------------------------
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
//
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901

#ifndef G4QIsotope_hh
#define G4QIsotope_hh

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <vector>

class G4QIsotope
{
protected:
  G4QIsotope(); // ***Singletone*** All natural elements are initialized in the Constructor

public:
  ~G4QIsotope(); // It's public for compilation purposes on Windows, user must not call it!

  // Create newElement with the Abundancy vector (User must delete elements of the vector)
  // -----------------------------------------------------------------------------------
  // Example of initialization of the new Element (not natural abanduncies of isotopes):
  // -----------------------------------------------------------------------------------
  //std::vector<std::pair<G4int,G4double> >*a= new std::vector<std::pair<G4int,G4double> >;
  // a->push_back(std::make_pair(n1,abundancy1));
  // a->push_back(std::make_pair(n2,abundancy2));
  // a->push_back(std::make_pair(n3,1.-abundancy1-abundancy2));
  ////Sum of abundancies must be 1 otherwise worning appears & if less LastAbu is increased
  // G4int Z=1;
  // G4int ind=1;
  //// ind>0, if =<0 then the InitElement member function returns the first free index !!!
  //// For ind>0 if theIndex "ind" already exists, returns the first free index        !!!
  // G4QIsotope::Get()->InitElement(Z, ind, a); // G4QIsotope class is a Singletone
  // std::for_each(a->begin(), a->end(), void operator()(std::pair<G4int,G4double> >* P)
  //                                                                         {delete P;});
  //// OR just use the following, which is faser:
  //// G4int nA=a->size();
  //// if(nA) for(G4int i=0; i<nA; i++) {delete a->operator[](i);}
  //
  G4int InitElement(G4int Z, G4int index, std::vector<std::pair<G4int,G4double>*>* abund);

  // The highest index defined for Element with Z (Index>0 correspondToUserDefinedElements)
  // -----------------------------------------------------------------------------------
  G4int GetLastIndex(G4int Z); // Returns the last defined index (if only natural: =0)

  // Indices can have differen numbers (not 1,2,3,...) & in different sequences (9,3,7,...)
  G4bool IsDefined(G4int Z, G4int Ind); // Returns true if defined, false, if not defined

  // A#ofNeutrons in Element with Z & UseDefIndex. Universal for Nat(index=0) & UserDefElem
  G4int GetNeutrons(G4int Z, G4int index=0);//If theElement doesn't exist, returns negative

  // #ofProtons in stable isotopes with fixed A=Z+N. Returns length and fils VectOfIsotopes
  // -----------------------------------------------------------------------------------
  // Example of printing of isotopes with A=152:
  // -------------------------------------------
  // G4int A=152;               // A can not be more than 269 
  // std::vector<G4int> isV(4); // At present A with nIso>4 are not known
  // G4int nIso= G4QIsotope::Get()->GetProtons(A, isV); // isV is cleaned up before filling
  // if(nIso)for(G4int i,i<nIso,i++)G4cout<<"I#"<<i<<"Z="<<isV[i]<<",N="<<A-isV[i]<<G4endl;
  G4int GetProtons(G4int A, std::vector<G4int>& isoV);

  // Get a pointer to the vector of pairs(N,CrosS), where N is used to calculate CrosS
  // -----------------------------------------------------------------------------------
  // Example of initialization of the Cross Section to randomize weighted isotopes:
  // -----------------------------------------------------------------------------------
  // std::vector<std::pair<G4int,G4double>*>* cs= G4QIsotope::Get()->GetCSVector(Z, index);
  // G4int nIs=cs->size; // A#Of Isotopes in the element
  // if(nIs) for(G4int i; i<nIs; i++)
  // {
  //   G4int N=cs->at(i)->first;   // A#Of neuterons in the isotope
  //   cs->at(i)->second = CalculateCrossSection(particle,Z,N)// Calc particle+A(Z,N) CrosS
  // }
  std::vector<std::pair<G4int,G4double>*>* GetCSVector(G4int Z, G4int index = 0);

  // Get the abundancy vector for calculation of mean cross sections
  std::vector<std::pair<G4int,G4double>*>* GetAbuVector(G4int Z, G4int index = 0);

  // Get the summed abundancy vector (e.g. for randomization by itself)
  std::vector<std::pair<G4int,G4double>*>* GetSumAVector(G4int Z, G4int index = 0);

  // Calculates the mean Cross Section for the initialized Element(ind=0 Nat,ind>0 UserDef)
  G4double GetMeanCrossSection(G4int Z, G4int index = 0); // IsoCS's must init, IfNotRet<0

  // Randomize A#OfNeutrons in the Isotope weighted by theAbubdancies and theCrossSections
  G4int GetCSNeutrons(G4int Z, G4int index = 0); // IsoCrosSections must init, IfNotRet<0

  static G4QIsotope* Get();         // Get a pointer to the Singletone G4QIsotope

private:
  G4int RandomizeNeutrons(G4int Z); // Gives a#of neutrons in the Random Isotope for the Z

private:
  // Initialized in the constructor
  static std::vector<std::vector<std::pair<G4int,G4double>*>*> natElements; //NaturalElem's
  static std::vector<std::vector<std::pair<G4int,G4double>*>*> natSumAbund; //NatElemSumA's
  static std::vector<std::vector<std::pair<G4int,G4double>*>*> natIsoCrosS; //CSOfNatElem's
  // It is initialized by user, but it is cleaned up in the destructor
  static std::vector<std::pair<G4int,std::vector<std::pair<G4int,G4double>*>*>*> newElems;
  static std::vector<std::pair<G4int,std::vector<std::pair<G4int,G4double>*>*>*> newSumAb;
  static std::vector<std::pair<G4int,std::vector<std::pair<G4int,G4double>*>*>*> newIsoCS;
};
#endif
