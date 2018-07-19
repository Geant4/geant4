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
// $Id: G4Evaporation.hh 102025 2016-12-16 14:43:41Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
//
// 14/02/2007 Alex Howard - added protection for negative probabilities in the sum 
// 27/07/2009 V.Ivanchenko - added Combined decay channels (default + GEM) 
// 11/05/2010 V.Ivanchenko - rewrited technical part do not "new" and "delete" 
//                           of small objects
// 22/04/2011 V.Ivanchenko - added maxZ and maxA for FermiBreakUp model
// 23/01/2012 V.Ivanchenko added pointer of G4VPhotonEvaporation to the constructor
// 06/05/2013 V.Ivanchenko added pointer to ion table
// 04/10/2014 D. Mancusi Moved theChannels and theChannelFactory to the base
//                       class, since they seem to be common to all classes
//                       derived from G4VEvaporation.
//

#ifndef G4Evaporation_h
#define G4Evaporation_h 1

#include "globals.hh"

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Fragment.hh"
#include "G4DeexPrecoParameters.hh"
#include <vector>

class G4VEvaporationFactory;
class G4NistManager;
class G4IonTable;
class G4VFermiBreakUp;
class G4UnstableFragmentBreakUp;

class G4Evaporation : public G4VEvaporation
{
public:

  explicit G4Evaporation(G4VEvaporationChannel* photoEvaporation = nullptr);
	 
  virtual ~G4Evaporation();

  virtual void InitialiseChannels() final;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products 
  virtual void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus) final;

  void SetDefaultChannel();
  void SetGEMChannel();
  void SetCombinedChannel();

private:

  void InitialiseChannelFactory();

  G4Evaporation(const G4Evaporation &right) = delete;
  const G4Evaporation & operator=(const G4Evaporation &right) = delete;
  G4bool operator==(const G4Evaporation &right) const = delete;
  G4bool operator!=(const G4Evaporation &right) const = delete;

  G4int    fVerbose;
  size_t   nChannels;
  G4double minExcitation;
  G4NistManager* nist;
  G4IonTable*    theTableOfIons;
  G4UnstableFragmentBreakUp* unstableBreakUp;
  G4bool isInitialised;

  G4DeexChannelType channelType;

  std::vector<G4double> probabilities;
};

#endif
