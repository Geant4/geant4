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
// $Id: G4DNACrossSectionDataSet.hh,v 1.1 2005-07-07 16:38:58 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Riccardo Capra <capra@ge.infn.it>
//
// History:
// -----------
// 30 Jun 2005  RC         Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef   G4DNACROSSSECTIONDATASET_HH
 #define  G4DNACROSSSECTIONDATASET_HH 1

 #include "G4ShellEMDataSet.hh"

 class G4DNACrossSectionDataSet : public G4ShellEMDataSet
 { 
  public:
                                                G4DNACrossSectionDataSet(G4VDataSetAlgorithm* argAlgorithm, G4double argUnitEnergies=MeV, G4double argUnitData=barn);
   virtual                                     ~G4DNACrossSectionDataSet();

   virtual G4bool                               LoadData(const G4String & argFileName);
   virtual G4bool                               SaveData(const G4String & argFileName) const;
   
  private:
   G4String                                     FullFileName(const G4String & argFileName) const;
 };
#endif /* G4DNACROSSSECTIONDATASET_HH */
