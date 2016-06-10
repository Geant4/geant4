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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4DPMJET2_5CrossSection.cc
/// \brief Implementation of the G4DPMJET2_5CrossSection class
//
// $Id: G4DPMJET2_5CrossSection.cc 81932 2014-06-06 15:39:45Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4DPMJET2_5CrossSection.cc
//
// Version:             0.A
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET


#include "G4DPMJET2_5CrossSection.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicException.hh"

#include "globals.hh"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "G4DynamicParticle.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
G4DPMJET2_5CrossSection::G4DPMJET2_5CrossSection ():
  upperLimit ( 1000.0 * TeV ), lowerLimit ( 5.0 * GeV ), maxA(240)
{
  theCrossSectionIndex.clear();
  Initialise();
  //
  //
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // This next bit is provisional, stating that this cross-section estimator
  // is applicable to hydrogen targets.  However, the cross-section will be
  // set to zero.
  //
  ATmin = 1;
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //
}
///////////////////////////////////////////////////////////////////////////////
//
G4DPMJET2_5CrossSection::~G4DPMJET2_5CrossSection ()
{
  //
  // Go through the list of cross-section fit parameters and delete the arrays.
  //
  G4cout << "G4DPMJET2_5CrossSection::~G4DPMJET2_5CrossSection" << G4endl;
  G4cout << "Size: " << theCrossSectionIndex.size() << G4endl;
  /*  
  if(theCrossSectionIndex.size() > 0) {

    G4DPMJET2_5CrossSectionIndex::iterator it;
    for (it=theCrossSectionIndex.begin(); it!=theCrossSectionIndex.end(); ++it)
      {
        G4DPMJET2_5CrossSectionParamSet *ptr = it->second;
        for (G4DPMJET2_5CrossSectionParamSet *ptr1=ptr; ptr1<ptr+maxA; ptr1++)
          { delete ptr1; }
      }
  }
  */
}
///////////////////////////////////////////////////////////////////////////////
//
G4bool 
G4DPMJET2_5CrossSection::IsIsoApplicable(const G4DynamicParticle* theProjectile, 
                                         G4int, G4int AT,
                                         const G4Element*, const G4Material*)
{
  G4bool result = false;
  G4int AP = theProjectile->GetDefinition()->GetBaryonNumber();
  if(AP >= 1) {
    G4double EPN = theProjectile->GetKineticEnergy()/G4double(AP);
    result = (EPN >= lowerLimit && EPN <= upperLimit &&
              AT  >= ATmin      && AT  <= ATmax &&
              AP  >= APmin      && AP  <= APmax);
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double 
G4DPMJET2_5CrossSection::GetIsoCrossSection(
    const G4DynamicParticle* theProjectile, 
    G4int ZZ, G4int AT, const G4Isotope*,
    const G4Element*, const G4Material*)
{
  //
  // Initialise the result.
  G4double result = 0.0;
  //
  //
  // Get details of the projectile and target (nucleon number, atomic number,
  // kinetic enery and energy/nucleon.
  //
  G4int AP     = theProjectile->GetDefinition()->GetBaryonNumber();
  G4double TP  = theProjectile->GetKineticEnergy();
  G4double EPN = TP /G4double(AP);

  if (AT < ATmin || AT > ATmax || AP < APmin || AP > APmax ||
      EPN < lowerLimit || EPN > upperLimit)
  {
    G4cout <<G4endl;
    G4cout <<"ERROR IN G4DPMJET2_5CrossSection::GetIsoZACrossSection" <<G4endl;
    G4cout <<"ATTEMPT TO USE CROSS-SECTION OUTSIDE OF RANGE"          <<G4endl;
    G4cout <<"NUCLEON NUMBER OF PROJECTILE = " <<AP                   <<G4endl;
    G4cout <<"NUCLEON NUMBER OF TARGET     = " <<AT                   <<G4endl;
    G4cout <<"ENERGY PER NUCLEON           = " <<EPN*MeV              <<G4endl;
    G4cout <<"ACCEPTABLE RANGE FOR AP      = " <<APmin
           <<" TO "                            <<APmax                <<G4endl;
    G4cout <<"ACCEPTABLE RANGE FOR AT      = " <<ATmin
           <<" TO "                            <<ATmax                <<G4endl;
    G4cout <<"ACCEPTABLE RANGE FOR ENERGY  = " <<lowerLimit
           <<" MeV/n TO "                      <<upperLimit
           <<" MeV/n" <<G4endl;
    G4cout <<G4endl;
    return result;
  }
  //
  //
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  // This next bit is provisional, stating that this cross-section hydrogen
  // targets is zero.
  //
  if ( AT == 1 ) return 0.0;
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //
  //
  //
  // Results are parameterised as a function of the natural logarithm of the
  // centre of mass energy of the projectile and target system.
  //
  G4double sigma = 0.0;
  G4double mT    = G4NucleiProperties::GetNuclearMass(AT, ZZ);
  G4double EP    = theProjectile->GetTotalEnergy();
  G4double mP    = theProjectile->GetDefinition()->GetPDGMass();
  G4double lnECM = std::log(std::sqrt(mP*mP + mT*mT + 2.0*mT*EP));
  G4DPMJET2_5CrossSectionIndex::iterator it = theCrossSectionIndex.find(AT);
  if (it != theCrossSectionIndex.end())
  {
    G4DPMJET2_5CrossSectionParamSet *ptr = (it->second) + AP;
    G4double cc0 = (*ptr)[0];
    G4double cc1 = (*ptr)[1];
    G4double cc2 = (*ptr)[2];
    sigma = cc0 + cc1*lnECM + cc2*lnECM*lnECM;
    sigma = sigma * millibarn;
    if (verboseLevel >= 2) {
      G4cout <<"***************************************************************"
             <<G4endl;
      G4cout <<"G4DPMJET2_5CrossSection::GetIsoCrossSection" <<G4endl;
      G4cout <<"PROJECTILE    = "
             <<theProjectile->GetDefinition()->GetParticleName() <<G4endl;
      G4cout <<"TARGET (A,Z)  = (" <<AT <<"," <<ZZ <<")" <<G4endl;
      G4cout <<"K. ENERGY/NUC = " <<EPN/MeV <<" MeV/n" <<G4endl;
      G4cout <<"CROSS SECTION = " <<sigma/millibarn <<" MILLIBARNS" <<G4endl;
      G4cout <<"***************************************************************"
             <<G4endl;
    }
  }
  else
  {
    G4cout <<G4endl;
    G4cout <<"ERROR IN G4DPMJET2_5CrossSection::GetIsoCrossSection" <<G4endl;
    G4cout <<"NO CROSS-SECTION FIT DATA LOADED FOR AT = " <<AT        <<G4endl;
    G4cout <<G4endl;
  }
  
  return sigma;
}
///////////////////////////////////////////////////////////////////////////////
//
void G4DPMJET2_5CrossSection::Initialise ()
{
  //verboseLevel = 2;
  //
  //
  // Determine first if the environment variable G4DPMJET2_5DATA is set.  If not
  // then ask for it to be set and call exception.
  //
  if ( !getenv("G4DPMJET2_5DATA") )
  {
    G4cout <<"ENVIRONMENT VARIABLE G4DPMJET2_5DATA NOT SET " <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
      "Please setenv G4DPMJET2_5DATA to point to the dpmjet2.5 data files.");
  }

  G4String filename = G4String(getenv("G4DPMJET2_5DATA")) + "/" +
    "GlauberCrossSections.dat";

  std::ifstream glauberXSFile(filename);
  if (glauberXSFile) {
    //
    //
    // Glaubercross-section file does exist, so read in maximum and minimum A
    // for target and projectile.
    //
    glauberXSFile >>APmin >>APmax >>ATmin >>ATmax;
    //
    //
    // Determine the list of targets based on the G4ElementList.  The list of
    // target nucleon numbers is stored as a ket to the map theCrossSectionIndex.
    // G4double[240][3] array objects are created to allow storage of the 
    // cross-section fit parameters.
    //
    const G4ElementTable *theElementTable = G4Element::GetElementTable();
    G4ElementTable::const_iterator it;
    for (it=theElementTable->begin(); it!=theElementTable->end(); it++)
      {
        G4int nIso = (*it)->GetNumberOfIsotopes();
        G4IsotopeVector* isoVector = (*it)->GetIsotopeVector();
        for (G4int i = 0; i < nIso; i++)
          {
            G4int AA = (*isoVector)[i]->GetN();
            if (theCrossSectionIndex.count(AA) == 0 && 
                AA >= ATmin && AA <= ATmax)
              {
//
//
// Whilst the use of std::map should eliminate duplication of keys, we need to
// know whether isotope's with the same nucleon number have been declared before
// creating the large arrays, hence the use of the "count" member function.
//
                G4DPMJET2_5CrossSectionParamSet *a =
                  new G4DPMJET2_5CrossSectionParamSet[maxA];
                theCrossSectionIndex.insert(
                  G4DPMJET2_5CrossSectionIndex::value_type(AA,a));
              }
          }
      }
    
    //
    //
    // Now proceed to read in the remainder of the GlauberCrossSection.dat file,
    // loading into theCrossSectionIndex any relevant fitting parameters to the
    // target nuclei.
    //
    char inputChars[80]={' '};
    G4String inputLine;
    while (-glauberXSFile.getline(inputChars, 80).eof() != EOF)
      {
        inputLine = inputChars;
        if (inputLine.length() != 0)
          {
            std::istringstream tmpStream(inputLine);
            G4int AP, AT;
            G4double cc0, cc1, cc2;
            tmpStream >>AP >>AT >>cc0 >>cc1 >>cc2;
            G4DPMJET2_5CrossSectionIndex::iterator IT = 
              theCrossSectionIndex.find(AT);
            if (IT != theCrossSectionIndex.end())
              {
                G4DPMJET2_5CrossSectionParamSet *ptr = (IT->second) + AP;
                *ptr = G4DPMJET2_5CrossSectionParamSet(cc0,cc1,cc2);
              }
          }
      }

    glauberXSFile.close();
    G4cout << "G4DPMJET2_5CrossSection::Initialise () done!" << G4endl;
  } else {
    G4cout <<"GlauberCrossSections.dat DOES NOT EXIST" <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
  "GlauberCrossSections.dat should be located in $G4DPMJET2_5DATA directory.");
  }
}
///////////////////////////////////////////////////////////////////////////////
//
void G4DPMJET2_5CrossSection::BuildPhysicsTable(const G4ParticleDefinition&)
{;}
///////////////////////////////////////////////////////////////////////////////
//
void G4DPMJET2_5CrossSection::DumpPhysicsTable(const G4ParticleDefinition 
  &theProjectile)
{
  const G4int AP    = G4lrint(theProjectile.GetBaryonNumber());
  G4cout <<G4endl;
  G4cout <<"G4DPMJET2_5CrossSection::DumpPhysicsTable" <<G4endl;
  G4cout <<"DUMPING CROSS-SECTION FITTING COEFFICIENTS FOR AP = "
         <<AP <<G4endl;
  G4cout <<G4endl;
  G4cout <<"   AT"
         <<"             c0"
         <<"             c1"
         <<"             c2"
         <<G4endl;
  G4DPMJET2_5CrossSectionIndex::iterator it;
  for (it=theCrossSectionIndex.begin(); it!=theCrossSectionIndex.end(); it++)
  {
    G4cout.unsetf(std::ios::scientific);
    G4cout.setf(std::ios::fixed|std::ios::right|std::ios::adjustfield);
    G4cout.precision(0);
    G4cout <<std::setw(5)  <<it->first;

    G4cout.unsetf(std::ios::fixed);
    G4cout.setf(std::ios::scientific|std::ios::right|std::ios::adjustfield);
    G4cout.precision(7);
    G4DPMJET2_5CrossSectionParamSet *ptr = (it->second) + AP;
    G4cout <<std::setw(15) <<(*ptr)[0]
           <<std::setw(15) <<(*ptr)[1]
           <<std::setw(15) <<(*ptr)[2]
           <<G4endl;
  }
  G4cout.setf(std::ios::fixed);
}
#endif
