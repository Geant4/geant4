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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPArbitaryTab.hh"
#include "G4ios.hh"

  G4double G4ParticleHPArbitaryTab::Sample(G4double anEnergy) 
  {
    G4int i;
    for(i=0;i<nDistFunc;i++)
    {
      if(anEnergy<theDistFunc[i].GetLabel()) break; // that is the energy we need
    }
    G4int low(0), high(0);
    if(i==nDistFunc) 
    {
      low = i-2;
      high = i-1;
    }
    else if(i==0)
    {
      if(nDistFunc==0)
      {
        G4cerr << "No distribution functions to sample "
             << "from in G4ParticleHPArbitaryTab::Sample"<<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "nDistFunc==0");
      } 
      else 
      {
        return theDistFunc[0].Sample();
      }
    }
    else
    {
      low = i-1;
      high = i;
    }
    //************************************************************************
    //EMendoza
    /*
      theBuffer.Merge(theManager.GetScheme(low), anEnergy, 
      theDistFunc+low, theDistFunc+high);
      return theBuffer.Sample();
    */
    //************************************************************************
    //New way to perform the 2D sampling:
    G4double elow=theDistFunc[low].GetLabel();
    G4double ehigh=theDistFunc[high].GetLabel();
    G4double rval=(anEnergy-elow)/(ehigh-elow);//rval is 0 for elow and 1 for ehigh
    G4double eoutlow=theLowThreshold[low]+rval*(theLowThreshold[high]-theLowThreshold[low]);
    G4double eouthigh=theHighThreshold[low]+rval*(theHighThreshold[high]-theHighThreshold[low]);
    G4double rand=G4UniformRand();
    G4double Eout_1=0,Eout_2=0;
    if(rval<rand){
      Eout_1=theDistFunc[low].Sample();
      Eout_2=eoutlow+(Eout_1-theLowThreshold[low])*(eouthigh-eoutlow)/(theHighThreshold[low]-theLowThreshold[low]);
    }
    else{
      Eout_1=theDistFunc[high].Sample();
      Eout_2=eoutlow+(Eout_1-theLowThreshold[high])*(eouthigh-eoutlow)/(theHighThreshold[high]-theLowThreshold[high]);
    }
    return Eout_2;

    //************************************************************************
  }
