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

#include "G4MSUHadData.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ElementVector.hh"

#include <fstream>
#include <strstream>
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//sim
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "G4HadFileFinder.hh"
//end sim

G4String G4MSUHadData::fMSUenergy = "MeV";
G4String G4MSUHadData::fMSUangle = "degree";
G4String G4MSUHadData::fMSUddXsc = "millibarn/sr/MeV";

G4MSUHadData::G4MSUHadData( G4String&  inputFn, G4HadFileSpec& filetowrite):
  G4VHadDataWriting()
{
  SetInputFileName(inputFn);
  SetEnergyUnit(fMSUenergy);
  SetAngleUnit(fMSUangle);
  SetDdXscUnit(fMSUddXsc);

  WriteDataFile(filetowrite);
}

G4MSUHadData::~G4MSUHadData()
{ 

}



///////////////////////////////////////////////////////////////
//
// Check dir/fn !!!

void G4MSUHadData::WriteDataFile( G4HadFileSpec& filetowrite) 
{

  G4String name(filetowrite.G4HDSFilename());
  G4String pathString(filetowrite.G4HDSFilepath());
  
  G4String dirFile = pathString + name;   
  // G4String dirFile = name;   

  char* doIhaverights = getenv("G4HADWRITTINGRIGHTS");

  int doI = strcmp(doIhaverights,"0");

  if ( ( !doIhaverights) || ( doI == 0 ) )
  {
      G4String excep = "G4HADWRITTINGRIGHTS environment variable is not set or equals 0. ";
      excep += "You don't have the permission to write in the G4HDS";
      //  G4Exception(excep);
      
  }
  G4cout << "Filling file in position " << dirFile << G4endl;

  FillDoubleDiffXSC(dirFile);

}

//////////////////////////////////////////////////////////////////
//
// Fill Tkin and xsc data for double differential xsc 
// per angle and per momentumC. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4MSUHadData::FillDoubleDiffXSC( G4String& fileName)
{
  G4String tmpString, tmpSs1, tmpSs2, tmpSs3, tmpSs4, tmpSs5, tmpSs6;
  size_t posValue, pos1, pos2;
  size_t position = 0;
  G4double Tkin, xsc, angle, energy, deltaE, deltaXsc;
  G4int angleNo, omegaNo;

  G4cout<<"G4VHadDataWriting::FillDifferentialXSC(fN),\n"<<" fN = "
        <<fileName<<G4endl;

  std::ifstream fin(fInputFileName);
  std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
    G4String excep = "G4HadDataSet - data file: " + fInputFileName + " not found";
    G4Exception(excep);
  }

  while( !fin.eof() )  // reading loop from a MSU file to vectors
  {
    getline(fin,tmpString);
    //  G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      fCommentVector->push_back(tmpString);
      //  G4cout<<tmpString<<G4endl;  // print comment lines
    }
    else if ( position == tmpString.find_first_of(fAnyEmptySpace) || 
               position == tmpString.find("\n")     || 
               position == tmpString.find("\r")        )
    {
      //  G4cout<<"tmpString is empty "<<G4endl;
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,position);
      //  G4cout<<"posValue = "<<posValue<<G4endl;
      if(posValue+1 >= tmpString.size() ) continue;      
      else if( posValue != 0) G4cout<<"Warning: file line starts with delims ?!"<<G4endl;
    }
    else if ( position == tmpString.find("energyUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      //  G4cout<<"energyUnit   "<<tmpSs1<<G4endl;
      SetEnergyUnit(tmpSs1);
    }
    else if ( position == tmpString.find("angleUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      // G4cout<<"angleUnit   "<<tmpSs1<<G4endl;
      SetAngleUnit(tmpSs1);
    }
    else if ( position == tmpString.find("ddXscUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      //  G4cout<<"ddXscUnit   "<<tmpSs1<<G4endl;
      SetDdXscUnit(tmpSs1);
    }
    else if ( position == tmpString.find("energyNo") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      fNo      = atoi(tmpSs1);
      fEnergyNoVector.push_back(fNo);
      //  G4cout<<"energyNo = fNo = "<<fNo<<G4endl;
    }
    else if ( position == tmpString.find("Tkin") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      Tkin = G4double(atof(tmpSs1));  
      fTkinVector->push_back(Tkin);
      //  G4cout<<"Tkin = "<<Tkin<<G4endl;
    }
    else if ( position == tmpString.find("angleNo") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      angleNo      = atoi(tmpSs1);
      fAngleNoVector.push_back(angleNo);
      //  G4cout<<"angleNo = "<<angleNo<<G4endl;
    }
    else if ( position == tmpString.find("angle") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      angle = G4double(atof(tmpSs1));  
      fAngleVector->push_back(angle);
      //  G4cout<<"angle = "<<angle<<G4endl;
    }
    else if ( position == tmpString.find("omegaNo") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      omegaNo      = atoi(tmpSs1);
      fOmegaNoVector.push_back(omegaNo);
      // G4cout<<"omegaNo = "<<omegaNo<<G4endl;
    }
    else // Tkin, ddXsc filling to class field vectors, errors=0
    { 

      pos1 = tmpString.find_first_of(fAnyEmptySpace);
      pos2 = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      //  pos3 = tmpString.find_first_of(fAnyEmptySpace,pos2);
      //  G4cout<<position<<"\t"<<pos1<<"\t"<<pos2<<"\t"<<pos3<<G4endl; 
  
      tmpSs1 = tmpString.substr(position,pos1);
      tmpSs2 = tmpString.substr(pos2);          // ,pos3-pos2);

      energy   = G4double(atof(tmpSs1));  
      xsc      = G4double(atof(tmpSs2));
      deltaE   = 0.;
      deltaXsc = 0.;        

      fMomentumCVector->push_back(energy);
      fXscVector      ->push_back(xsc);       
      fDeltaMomCVector->push_back(deltaE);
      fDeltaXscVector ->push_back(deltaXsc);       

      // G4cout<<energy<<"\t"<<deltaE<<"\t"<<xsc<<"\t"<<deltaXsc<<G4endl;
      //  G4cout<<pos3-pos2<<"\t";
    }
  }
  fin.close();

  G4cout<<"reading closed; writing up"<<G4endl;
  // VG: so the data were stored in vectors and now we are ready fileout 
  // them in G4HDS style. We introduce some integers to struct the data

  size_t iComment, iEnergy, jAngle, kOmega; 
  size_t k = 0, j = 0; 

  struct stat s;
  G4int hfile=stat(fileName,&s);
  if(hfile !=-1)    G4cout<<"I found you"<< G4endl;
      
  
  std::ofstream fout(fileName, std::ios::out | std::ios::app );
  fout.setf( std::ios::scientific, std::ios::floatfield );

  // start writing from memory vectors to G4HDS

  for( iComment = 0; iComment < fCommentVector->size(); ++iComment)
  {
    fout<<(*fCommentVector)[iComment]<<G4endl;
    G4cout<<(*fCommentVector)[iComment]<<G4endl;
  }
  fout<<G4endl<<"energyNo"<<"\t"<<fNo<<G4endl;
  G4cout<<G4endl<<"energyNo"<<"\t"<<fNo<<G4endl;
 

  for(iEnergy = 0; iEnergy < fNo; ++iEnergy)
  {
    fout<<G4endl<<"Tkin"<<"\t"<<(*fTkinVector)[iEnergy]<<G4endl;
    G4cout<<G4endl<<"Tkin"<<"\t"<<(*fTkinVector)[iEnergy]<<G4endl;

    fout<<G4endl<<"energyUnit"<<"\t"<<fEnergyUnit<<G4endl;
    fout<<G4endl<<"angleUnit"<<"\t"<<fAngleUnit<<G4endl;
    fout<<G4endl<<"ddXscUnit"<<"\t"<<fDdXscUnit<<G4endl;
    G4cout<<G4endl<<"energyUnit"<<"\t"<<fEnergyUnit<<G4endl;
    G4cout<<G4endl<<"angleUnit"<<"\t"<<fAngleUnit<<G4endl;
    G4cout<<G4endl<<"ddXscUnit"<<"\t"<<fDdXscUnit<<G4endl;
    
    fout<<G4endl<<"angleNo"<<"\t"<<fAngleNoVector[iEnergy]<<G4endl<<G4endl;
    G4cout<<G4endl<<"angleNo"<<"\t"<<fAngleNoVector[iEnergy]<<G4endl<<G4endl;

    for(jAngle = 0; jAngle < fAngleNoVector[iEnergy]; ++jAngle)
    {
      fout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j]<<G4endl<<G4endl;
      G4cout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j]<<G4endl<<G4endl;

      fout<<G4endl<<"omegaNo"<<"\t"<<fOmegaNoVector[j]<<G4endl<<G4endl;
      G4cout<<G4endl<<"omegaNo"<<"\t"<<fOmegaNoVector[j]<<G4endl<<G4endl;

      for(kOmega = 0; kOmega < fOmegaNoVector[j]; ++kOmega)
      {
        fout<<(*fMomentumCVector)[k]<<"\t"
                    <<(*fDeltaMomCVector)[k]<<"\t"
                    <<(*fXscVector)[k]<<"\t"
	            <<(*fDeltaXscVector)[k]<<G4endl;
        G4cout<<(*fMomentumCVector)[k]<<"\t"
                      <<(*fDeltaMomCVector)[k]<<"\t"
                      <<(*fXscVector)[k]<<"\t"
	              <<(*fDeltaXscVector)[k]<<G4endl;
        k++;
      }
      j++;
    } 
  }
  
  // sim Interpolation scheme 
  
  fout << G4endl << "# Space left for data interpolation scheme declaration"
       << G4endl << G4endl;
  G4cout << G4endl << "# Space left for data interpolation scheme declaration"
       << G4endl << G4endl;

  fout<<G4endl<<"# End data file" <<G4endl<<G4endl;
  G4cout<<G4endl<<"# End data file" <<G4endl<<G4endl;
}










