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

#include "G4HadDataReading.hh"
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

static const G4String G4HadDataReading::fAnyNumber = "1234567890.+-";
static const G4String G4HadDataReading::fAnyEmptySpace = " \t";
static const G4String G4HadDataReading::fAnyHidden = " \t\n\r";


G4HadDataReading::G4HadDataReading()
{
  fEnergyUnit     = MeV; 
  fAngleUnit      = degree; 
  fXscUnit        = millibarn;
  fXscPerAngleUnit = millibarn/sr; 
  fXscPerMomCUnit = millibarn/GeV;
  fDdXscUnit      = millibarn/sr/MeV;

  fTkinBin = -1.;
  fNo = 0;
  fTkinVector          = new G4DataVector();  
  fTkinBinVector       = new G4DataVector();  
  fXscVector           = new G4DataVector();  
  fDeltaXscVector      = new G4DataVector();  
  fMultiplicityVector  = new G4DataVector();  
  fMomentumCVector     = new G4DataVector;
  fDeltaMomCVector     = new G4DataVector;
  fMomentumCBinVector  = new G4DataVector;
  fAngleVector         = new G4DataVector();
  fAngleBinVector      = new G4DataVector();

  fAngleTable             = new G4PhysicsTable();
  fAngleDdTable             = new std::vector<G4DataVector*>;
  fMomentumCTable         = new G4PhysicsTable();
  fDoubleDiffXscBank      = new std::vector<G4PhysicsTable*>;
  fDoubleDiffXscErrorBank = new std::vector<G4PhysicsTable*>;
  fCommentVector          = new std::vector<G4String>;
}



G4HadDataReading::~G4HadDataReading()
{ 
  delete fTkinVector;
  delete fTkinBinVector;

  delete fXscVector;
  delete fDeltaXscVector;
  delete fMultiplicityVector;

  delete fMomentumCVector;
  delete fDeltaMomCVector;
  delete fMomentumCBinVector;

  delete fAngleVector;
  delete fAngleBinVector;

  fAngleTable->clearAndDestroy();
  fMomentumCTable->clearAndDestroy();

  delete fDoubleDiffXscBank;
}


//////////////////////////////////////////////////////////////////
//
// Load Tkin and xsc data for total, elastic and inelastic xsc. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4HadDataReading::LoadIntegralXSC( G4String& fileName)
{
  G4cout<<"G4HadDataReading::LoadIntegralXSC(fN),\n"<<" fN = "<<fileName<<G4endl;

  std::ifstream fin(fileName);
  std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
      G4String excep = "G4HadDataSet - data file: " + fileName + " not found";
      G4Exception(excep);
  }
  G4String tmpString;
  G4String tmpSs1;
  G4String tmpSs2;
  size_t position = 0;
  size_t posValue;
  G4double Tkin, xsc;

  while( !fin.eof() )
  {
    getline(fin,tmpString);
    // G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      G4cout<<tmpString<<G4endl;  // print comment lines
    }
    else if ( position == tmpString.find("bin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fTkinBin = atof(tmpSs1)*GeV;
      G4cout<<"fTkinBin = "<<fTkinBin/GeV<<G4endl;
    }
    else if ( position == tmpString.find("no") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fNo = atoi(tmpSs1);
      G4cout<<"fNo = "<<fNo<<G4endl;
    }
    else if ( tmpString.empty() ) // just skip empty lines
    {
      //     G4cout<<"tmpString is empty "<<G4endl;
      continue;
    }
    else // Tkin, xsc loading to class field vectors
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(0,posValue);
      tmpSs2   = tmpString.substr(posValue+1);

      Tkin     = atof(tmpSs1)*GeV;
      xsc      = atof(tmpSs2)*millibarn;

      fTkinVector->push_back(Tkin);
      fXscVector->push_back(xsc);

      G4cout<<"Tkin = "<<Tkin/GeV<<"\t xsc = "<<xsc/millibarn<<G4endl;
    }
  }
  if ( fNo != (G4int)fTkinVector->size() )
  {
    G4cout<<" fNo != fTkinVector->size(). Change no value in file"<<G4endl;
    G4cout<<" now fNo = fTkinVector->size()"<<G4endl;
    fNo = (G4int)fTkinVector->size();
  }
  for( size_t i = 0 ; i < fTkinVector->size(); ++i)
  {
    G4cout<<(*fTkinVector)[i]/GeV<<"\t"<<(*fXscVector)[i]/millibarn<<G4endl;
  }
  fin.close();  
}



//////////////////////////////////////////////////////////////////
//
// Load Tkin and multiplicity(Tkin). 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4HadDataReading::LoadMultiplicity( G4String& fileName)
{
  G4cout<<"G4HadDataReading::LoadIntegralXSC(fN),\n"<<" fN = "<<fileName<<G4endl;

  std::ifstream fin(fileName);
  std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
      G4String excep = "G4HadDataSet - data file: " + fileName + " not found";
      G4Exception(excep);
  }
  G4String tmpString;
  G4String tmpSs1;
  G4String tmpSs2;
  size_t position = 0;
  size_t posValue;
  G4double Tkin, mult;

  while( !fin.eof() )
  {
    getline(fin,tmpString);
    // G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      G4cout<<tmpString<<G4endl;  // print comment lines
    }
    else if ( position == tmpString.find("bin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fTkinBin = atof(tmpSs1)*GeV;
      G4cout<<"fTkinBin = "<<fTkinBin/GeV<<G4endl;
    }
    else if ( position == tmpString.find("no") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fNo = atoi(tmpSs1);
      G4cout<<"fNo = "<<fNo<<G4endl;
    }
    else if ( tmpString.empty() ) // just skip empty lines
    {
      //     G4cout<<"tmpString is empty "<<G4endl;
      continue;
    }
    else // Tkin, xsc loading to class field vectors
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(0,posValue);
      tmpSs2   = tmpString.substr(posValue+1);

      Tkin     = atof(tmpSs1)*GeV;
      mult      = atof(tmpSs2);

      fTkinVector->push_back(Tkin);
      fMultiplicityVector->push_back(mult);

      G4cout<<"Tkin = "<<Tkin/GeV<<"\t multiplicity = "<<mult<<G4endl;
    }
  }
  if ( fNo != (G4int)fTkinVector->size() )
  {
    G4cout<<" fNo != fTkinVector->size(). Change no value in file"<<G4endl;
    G4cout<<" now fNo = fTkinVector->size()"<<G4endl;
    fNo = (G4int)fTkinVector->size();
  }
  for( size_t i = 0 ; i < fTkinVector->size(); ++i)
  {
    G4cout<<(*fTkinVector)[i]/GeV<<"\t"<<(*fMultiplicityVector)[i]<<G4endl;
  }
  fin.close();  
}


//////////////////////////////////////////////////////////////////
//
// Load Tkin and xsc data for differential xsc per angle or per momentumC. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4HadDataReading::LoadDifferentialXSC( G4String& fileName, G4bool momORangle)
{
  G4cout<<"G4HadDataReading::LoadDifferentialXSC(fN),\n"<<" fN = "
        <<fileName<<G4endl;

  std::ifstream fin(fileName);
  std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
    G4String excep = "G4HadDataSet - data file: " + fileName + " not found";
    G4Exception(excep);
  }
  G4String tmpString, tmpSs1, tmpSs2;
  size_t posValue, position = 0;
  G4double Tkin, variable, xsc;
  G4int iVector;
  G4PhysicsFreeVector* dataVector = NULL;

  while( !fin.eof() )
  {
    getline(fin,tmpString);
    // G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      G4cout<<tmpString<<G4endl;  // print comment lines
    }
    else if ( position == tmpString.find("bin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fTkinBin = atof(tmpSs1)*GeV;
      G4cout<<"fTkinBin = "<<fTkinBin/GeV<<G4endl;
    }
    else if ( position == tmpString.find("no") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fNo = atoi(tmpSs1);
      iVector = 0;
      dataVector = new G4PhysicsFreeVector(fNo);
      G4cout<<"fNo = "<<fNo<<G4endl;
    }
    else if ( position == tmpString.find("Tkin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      Tkin = atof(tmpSs1)*GeV;
      fTkinVector->push_back(Tkin);
      G4cout<<"Tkin = "<<Tkin/GeV<<G4endl;
    }
    else if ( tmpString.empty() ) // just skip empty lines
    {
      //     G4cout<<"tmpString is empty "<<G4endl;
      continue;
    }
    else // d(xsc)/d(variable) at a given Tkin are loaded to class field vectors
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(0,posValue);
      tmpSs2   = tmpString.substr(posValue+1);

      if(momORangle)  // per momentumC
      {
        variable     = atof(tmpSs1)*GeV;
        xsc      = atof(tmpSs2)*millibarn/GeV;
        dataVector->PutValue(iVector,variable,xsc);
        iVector++;
        if(iVector == fNo)
	{
          fMomentumCTable->push_back(dataVector);
          delete dataVector; 
	} 
      }
      else   // per angle
      {
        variable     = atof(tmpSs1)*degree;        // *rad, *degree ? 
        xsc      = atof(tmpSs2)*millibarn/degree;
        dataVector->PutValue(iVector,variable,xsc);
        iVector++;
        if(iVector == fNo)
	{
	  fAngleTable->push_back(dataVector);
          delete dataVector; 
	}       
      }
      G4cout<<"variable = "<<variable<<"\t xsc = "<<xsc/millibarn<<G4endl;
    }
  }
  fin.close();
 
  if(momORangle)  // per momentumC
  {
    for(  size_t i = 0; i < fMomentumCTable->length(); ++i )
    {
      for(size_t k = 0; k < (*fMomentumCTable)(i)->GetVectorLength(); ++k)
      {
        G4cout<<(*fMomentumCTable)(i)->GetLowEdgeEnergy(k)<<"\t"
	      <<(*(*fMomentumCTable)(i))(k)<<G4endl;
      }
    } 
  }
  else  // per angle
  {
    for(  size_t i = 0; i < fAngleTable->length(); ++i )
    {
      for(size_t k = 0; k < (*fAngleTable)(i)->GetVectorLength(); ++k)
      {
	  G4cout<<(*fAngleTable)(i)->GetLowEdgeEnergy(k)<<"\t"
	    <<(*(*fAngleTable)(i))(k)<<G4endl;
      }
    } 
  }


}

//////////////////////////////////////////////////////////////////
//
// Load Tkin and xsc data for double differential xsc 
// per angle and per momentumC. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4HadDataReading::LoadDoubleDiffXSC( G4String& fileName)
{
  G4String tmpString, tmpSs1,tmpSs2,tmpSs3, tmpSs4;
  size_t posValue, pos1, pos2, pos3, pos4, pos5, pos6, pos7; 
  size_t position = 0;
  G4double Tkin, xsc, angle, energy, deltaE, deltaXsc;
  G4int angleNo, omegaNo=0;
  G4PhysicsFreeVector* dataVector = NULL;
  G4PhysicsFreeVector* errorVector = NULL;
  G4PhysicsTable*  dataTable = NULL;
  G4PhysicsTable*  errorTable = NULL;
  G4DataVector*    angleVector = NULL;

  G4cout<<"G4HadDataReading::LoadDifferentialXSC(fN),\n"<<" fN = "
        <<fileName<<G4endl;

  std::ifstream fin(fileName);
  std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
    G4String excep = "G4HadDataSet - data file: " + fileName + " not found";
    G4Exception(excep);
  }
  fNo=0;
  while( !fin.eof() )
  {
    getline(fin,tmpString);
    //  G4cout<<tmpString<<G4endl;
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      fCommentVector->push_back(tmpString);
      //  G4cout<<tmpString<<G4endl;  // print comment lines
    }
    else if ( tmpString.empty() ) continue;
    //  else if ( position == tmpString.find_first_of(fAnyHidden)) continue;

    /*    else if ( position == tmpString.find_first_of(fAnyEmptySpace) || 
               position == tmpString.find("\n")     || 
               position == tmpString.find("\r")        )
    {
      //  G4cout<<"tmpString is empty "<<G4endl;
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,position);
      //  G4cout<<"posValue = "<<posValue<<G4endl;
      if(posValue+1 >= tmpString.size() ) continue;      
      else if( posValue != 0) G4cout<<"Warning: file line starts with delims ?!"<<G4endl;
    }
    */
    else if ( position == tmpString.find("energyUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      //  G4cout<<"energyUnit   "<<tmpSs1<<G4endl;
      SetEnergyUnit(tmpSs1);
      fEnergyUnitVector.push_back(fEnergyUnit);
    }
    else if ( position == tmpString.find("angleUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      // G4cout<<"angleUnit   "<<tmpSs1<<G4endl;
      SetAngleUnit(tmpSs1);
      fAngleUnitVector.push_back(fAngleUnit);
    }
    else if ( position == tmpString.find("ddXscUnit") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      //  G4cout<<"ddXscUnit   "<<tmpSs1<<G4endl;
      SetDdXscUnit(tmpSs1);
      fDdXscUnitVector.push_back(fDdXscUnit);
    }
    else if ( position == tmpString.find("energyNo") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      fNo      += atoi(tmpSs1);
      fEnergyNoVector.push_back(fNo);
      //  G4cout<<"energyNo = fNo = "<<fNo<<G4endl;
    }
    else if ( position == tmpString.find("Tkin") )
    {
      pos1     = tmpString.find_first_of(fAnyEmptySpace);
      posValue = tmpString.find_first_not_of(fAnyEmptySpace,pos1);
      tmpSs1   = tmpString.substr(posValue);
      Tkin = G4double(atof(tmpSs1))*fEnergyUnit;  
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
      angle = G4double(atof(tmpSs1))*fAngleUnit;  
      fAngleVector->push_back(angle);
      //  fOmegaNoVector.push_back(omegaNo);
      omegaNo = 0;
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
 
    /*
    else if ( position == tmpString.find("enddata") )
    {
      fOmegaNoVector.push_back(omegaNo);
      // G4cout<<"last omegaNo = "<<omegaNo<<G4endl;
    }
    */

    else // Tkin, ddXsc filling to class field vectors, errors=0
    { 
      // omegaNo++;

      pos1 = tmpString.find_first_of(fAnyNumber);
      pos2 = tmpString.find_first_of(fAnyEmptySpace,pos1);
      pos3 = tmpString.find_first_of(fAnyNumber,pos2);
      pos4 = tmpString.find_first_of(fAnyEmptySpace,pos3);
      pos5 = tmpString.find_first_of(fAnyNumber,pos4);
      pos6 = tmpString.find_first_of(fAnyEmptySpace,pos5);
      pos7 = tmpString.find_first_of(fAnyNumber,pos6);

      // G4cout<<pos1<<"  "<<pos2<<"  "<<pos3<<"  "
      //       <<pos4<<"  "<<pos5<<"  "<<pos6<<"  "
      //       <<pos7<<G4endl; 
  
      tmpSs1 = tmpString.substr(pos1,pos2-pos1);
      tmpSs2 = tmpString.substr(pos3,pos4-pos3);
      tmpSs3 = tmpString.substr(pos5,pos6-pos5);
      tmpSs4 = tmpString.substr(pos7);
         

      energy   = G4double(atof(tmpSs1))*fEnergyUnit;  
      deltaE   = G4double(atof(tmpSs2))*fEnergyUnit;
      xsc      = G4double(atof(tmpSs3))*fDdXscUnit;
      deltaXsc = G4double(atof(tmpSs4))*fDdXscUnit;        

      fMomentumCVector->push_back(energy);
      fXscVector      ->push_back(xsc);       
      fDeltaMomCVector->push_back(deltaE);
      fDeltaXscVector ->push_back(deltaXsc);       

      // G4cout<<energy/fEnergyUnit<<"\t"<<deltaE/fEnergyUnit<<"\t"
      //       <<xsc/fDdXscUnit<<"\t"<<deltaXsc/fDdXscUnit<<G4endl;
   
    }
    
  }
  fin.close();

  size_t iComment, iEnergy, jAngle, kOmega; 
  size_t k = 0, j = 0; 

  for( iComment = 0; iComment < fCommentVector->size(); ++iComment)
  {
    // G4cout<<(*fCommentVector)[iComment]<<G4endl;
  }
  G4cout<<G4endl<<"total energyNo in file "<<"\t"<<fNo<<G4endl;
 

  for(iEnergy = 0; iEnergy < fNo; ++iEnergy)
  {
    //  G4cout<<G4endl<<"Tkin"<<"\t"<<(*fTkinVector)[iEnergy]/fEnergyUnit<<G4endl;

    //  G4cout<<G4endl<<"energyUnit"<<"\t"<<fEnergyUnitVector[iEnergy]<<G4endl;
    // G4cout<<G4endl<<"angleUnit"<<"\t"<<fAngleUnitVector[iEnergy]<<G4endl;
    // G4cout<<G4endl<<"ddXscUnit"<<"\t"<<fDdXscUnitVector[iEnergy]<<G4endl;
    
    //  G4cout<<G4endl<<"angleNo"<<"\t"<<fAngleNoVector[iEnergy]<<G4endl<<G4endl;

    dataTable   = new G4PhysicsTable(fAngleNoVector[iEnergy]);
    errorTable  = new G4PhysicsTable(fAngleNoVector[iEnergy]);
    angleVector = new G4DataVector(); 

    for(jAngle = 0; jAngle < fAngleNoVector[iEnergy]; ++jAngle)
    {
      //  G4cout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j]/fAngleUnitVector[iEnergy]
      //        <<G4endl<<G4endl;
      //  G4cout<<G4endl<<"omegaNo"<<"\t"<<fOmegaNoVector[j]<<G4endl<<G4endl;

      dataVector  =  new G4PhysicsFreeVector(fOmegaNoVector[j]);
      errorVector =  new G4PhysicsFreeVector(fOmegaNoVector[j]);

      for(kOmega = 0; kOmega < fOmegaNoVector[j]; ++kOmega)
      {
	//   G4cout<<(*fMomentumCVector)[k]/fEnergyUnitVector[iEnergy]<<"\t"
        //      <<(*fDeltaMomCVector)[k]/fEnergyUnitVector[iEnergy]<<"\t"
	//       <<(*fXscVector)[k]/fDdXscUnitVector[iEnergy]<<"\t"
	//      <<(*fDeltaXscVector)[k]/fDdXscUnitVector[iEnergy]<<G4endl;
       
        dataVector ->PutValue(kOmega,(*fMomentumCVector)[k],(*fXscVector)[k]);
        errorVector->PutValue(kOmega,(*fDeltaMomCVector)[k],(*fDeltaXscVector)[k]);

        k++;
      }
      dataTable ->push_back(dataVector);
      errorTable->push_back(errorVector);
      angleVector->push_back((*fAngleVector)[j]);
      j++;
    }
    fDoubleDiffXscBank->push_back(dataTable); 
    fDoubleDiffXscErrorBank->push_back(errorTable); 
    fAngleDdTable->push_back(angleVector);
  }


  /*
  G4cout<<G4endl<<"Output of fDoubleDiffXscBank data"<<G4endl;
  G4cout<<"fDoubleDiffXscBank->size() = "<<fDoubleDiffXscBank->size()<<G4endl<<G4endl;
  jAngle=0; 
  for(size_t i = 0; i < fDoubleDiffXscBank->size(); ++i)
  {
    G4cout<<G4endl<<"(*fDoubleDiffXscBank)["<<i<<"]->size() = "
        <<(*fDoubleDiffXscBank)[i]->size()<<G4endl;

    for(  size_t j = 0; j < (*fDoubleDiffXscBank)[i]->size(); ++j )
    {
      G4cout<<G4endl<<"(*(*fDoubleDiffXscBank)["<<i<<"])("<<j<<")->GetVectorLength() = "
        <<(*(*fDoubleDiffXscBank)[i])(j)->GetVectorLength()<<G4endl;
      G4cout<<G4endl<<"Tkin"<<"\t"<<"angle"<<"\t"
            <<"omega"<<"\t"<<"ddXsc"<<G4endl<<G4endl;

      for(size_t k = 0; k < (*(*fDoubleDiffXscBank)[i])(j)->GetVectorLength(); ++k)
      {
        G4cout<<(*fTkinVector)[i]/fEnergyUnitVector[i]<<"\t"
              <<(*fAngleVector)[jAngle]/fAngleUnitVector[i]<<"\t"
              <<(*(*fDoubleDiffXscBank)[i])(j)->GetLowEdgeEnergy(k)/fEnergyUnitVector[i]
              <<"\t"
	      <<(*(*(*fDoubleDiffXscBank)[i])(j))(k)/fDdXscUnitVector[i]<<G4endl;
      }
      jAngle++;
    } 
  }
  */

}

//////////////////////////////////////////////////////////
//
// Check exact unit string !!!

void G4HadDataReading::SetEnergyUnit(G4String& unitString)
{
  if      ( unitString == "MeV" ) fEnergyUnit = MeV;
  else if ( unitString == "keV" ) fEnergyUnit = keV;
  else if ( unitString == "eV" )  fEnergyUnit = eV;
  else if ( unitString == "TeV" ) fEnergyUnit = TeV;
  else if ( unitString == "PeV" ) fEnergyUnit = PeV;
}

//////////////////////////////////////////////////////////
//
// Check units first. It is sensitive to exact unit string !

void G4HadDataReading::SetAngleUnit(G4String& unitString)
{
  if      ( unitString == "radian" ||
            unitString == "rad"  )  fAngleUnit = rad;
  else if ( unitString == "mrad" ) fAngleUnit = mrad;
  else if ( unitString == "sr" )   fAngleUnit = sr;
  else if ( unitString == "deg" ||
            unitString == "degree")  fAngleUnit = degree;
}

//////////////////////////////////////////////////////////

void G4HadDataReading::SetXscUnit(G4String& unitString)
{
  if      ( unitString == "barn" )      fXscUnit = barn;
  else if ( unitString == "microbarn" ) fXscUnit = microbarn;
  else if ( unitString == "nanobarn" )  fXscUnit = nanobarn;
  else if ( unitString == "picobarn" )  fXscUnit = picobarn;
  else if ( unitString == "millibarn" ||
            unitString == "mb")  fXscUnit = millibarn;
}

//////////////////////////////////////////////////////////
//
// Check exact unit string !!!

void G4HadDataReading::SetDdXscUnit(G4String& unitString)
{
  if      ( unitString == "barn/MeV/sr" ||
            unitString == "barn/sr/MeV"  )       fDdXscUnit = barn/MeV/sr;

  else if ( unitString == "microbarn/MeV/sr" ||
            unitString == "microbarn/sr/MeV"   ) fDdXscUnit = microbarn/MeV/sr;

  else if ( unitString == "nanobarn/MeV/sr" ||
            unitString == "nanobarn/sr/MeV"    ) fDdXscUnit = nanobarn/MeV/sr;

  else if ( unitString == "picobarn/MeV/sr" ||
            unitString == "picobarn/sr/MeV"    ) fDdXscUnit = picobarn/MeV/sr;

  else if ( unitString == "millibarn/MeV/sr" ||
            unitString == "millibarn/sr/MeV" ||
            unitString == "mb/sr/MeV"        ||
            unitString == "mb/MeV/sr"          ) fDdXscUnit = millibarn/MeV/sr;
}



//
//
/////////////////////////////////////////////////////////////



