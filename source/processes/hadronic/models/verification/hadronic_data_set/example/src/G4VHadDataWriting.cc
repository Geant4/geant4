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

#include "G4VHadDataWriting.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4ElementVector.hh"

#include "g4std/fstream"
#include "g4std/strstream"
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//sim
#include <sys/stat.h>
//end sim

G4VHadDataWriting::G4VHadDataWriting()
{
  fEnergyUnit     = "GeV"; 
  fAngleUnit      = "degree"; 
  fXscUnit        = "millibarn";
  fXscPeAngleUnit = "millibarn/degree"; 
  fXscPerMomCUnit = "millibarn/GeV";
  fDdXscUnit      = "millibarn/degree/GeV";

  fTkinBin = 1.;
  fNo = 1;
  fTkinVector          = new G4DataVector();  
  fTkinBinVector       = new G4DataVector();  
  fXscVector           = new G4DataVector();  
  fMultiplicityVector  = new G4DataVector();  
  fMomentumCVector     = new G4DataVector;
  fMomentumCBinVector  = new G4DataVector;
  fAngleVector         = new G4DataVector();
  fAngleBinVector      = new G4DataVector();

  fAngleTable        = new G4PhysicsTable();
  fMomentumCTable    = new G4PhysicsTable();
  fDoubleDiffXscBank = new G4std::vector<G4PhysicsTable*>;
  fCommentVector     = new G4std::vector<G4String>;
}



G4VHadDataWriting::~G4VHadDataWriting()
{ 
  delete fTkinVector;
  delete fTkinBinVector;

  delete fXscVector;
  delete fMultiplicityVector;

  delete fMomentumCVector;
  delete fMomentumCBinVector;

  delete fAngleVector;
  delete fAngleBinVector;

  fAngleTable->clearAndDestroy();
  fMomentumCTable->clearAndDestroy();

  delete fDoubleDiffXscBank;
  delete fCommentVector;
}


//////////////////////////////////////////////////////////////////
//
// Fill Tkin and xsc data for total, elastic and inelastic xsc. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4VHadDataWriting::FillIntegralXSC( G4String& fileName)
{
  G4cout<<"G4VHadDataWriting::FillIntegralXSC(fN),\n"<<" fN = "<<fileName<<G4endl;

  G4std::ifstream fin(fInputFileName);
  G4std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
    G4String excep = "G4HadDataSet - data file: " + fInputFileName + " not found";
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
    else // Tkin, xsc Filling to class field vectors
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
// Fill Tkin and multiplicity(Tkin). 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4VHadDataWriting::FillMultiplicity( G4String& fileName)
{
  G4cout<<"G4VHadDataWriting::FillIntegralXSC(fN),\n"<<" fN = "<<fileName<<G4endl;

  G4std::ifstream fin(fileName);
  G4std::filebuf* finbuf = fin.rdbuf();
  
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
    else // Tkin, xsc Filling to class field vectors
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
// Fill Tkin and xsc data for differential xsc per angle or per momentumC. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4VHadDataWriting::FillDifferentialXSC( G4String& fileName, G4bool momORangle)
{
  G4cout<<"G4VHadDataWriting::FillDifferentialXSC(fN),\n"<<" fN = "
        <<fileName<<G4endl;

  G4std::ifstream fin(fileName);
  G4std::filebuf* finbuf = fin.rdbuf();
  
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
    else // d(xsc)/d(variable) at a given Tkin are Filled to class field vectors
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
// Fill Tkin and xsc data for double differential xsc 
// per angle and per momentumC. 
// Lines start with comments (#,//,*), "bin", "no" or numbers (Tkin, "\t", xsc).
// Lines with comments are printed.
// Lines with "bin", "no" followed by "\t" and corresponding numbers.
// Empty lines are skipped 

void G4VHadDataWriting::FillDoubleDiffXSC( G4String& fileName)
{
  G4cout<<"G4VHadDataWriting::FillDifferentialXSC(fN),\n"<<" fN = "
        <<fileName<<G4endl;

  G4std::ifstream fin(fInputFileName);
  G4std::filebuf* finbuf = fin.rdbuf();
  
  if ( !( finbuf->is_open() ) )
  {
    G4String excep = "G4HadDataSet - data file: " + fInputFileName + " not found";
    G4Exception(excep);
  }
  G4String tmpString, tmpSs1, tmpSs2, tmpSs3, tmpSs4, tmpSs5, tmpSs6;
  size_t posValue, position = 0, pos1, pos2, pos3, pos4,pos5;
  G4double Tkin, variable, xsc, angleBin, angle, energy, deltaE, deltaXsc;
  G4int kTable, jVector=0, angleNo;
  G4PhysicsFreeVector* dataVector = NULL;
  G4PhysicsTable*  dataTable = NULL;

  while( !fin.eof() )
  {
    getline(fin,tmpString);
    // G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      fCommentVector->push_back(tmpString);
      G4cout<<tmpString<<G4endl;  // print comment lines
    }
    /* ************************************************************

    else if ( position == tmpString.find("angleUnit") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      SetAngleUnit(tmpSs1);
    }
    else if ( position == tmpString.find("energyUnit") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      SetEnergyUnit(tmpSs1);
    }
    else if ( position == tmpString.find("ddXscUnit") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      SetDdXscUnit(tmpSs1);
    }
    else if ( position == tmpString.find("bin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fTkinBin = atof(tmpSs1)*fEnergyUnit;
      G4cout<<"fTkinBin = "<<fTkinBin/GeV<<G4endl;
    }
    // Treatment of different angles at a given energy

    else if ( position == tmpString.find("angleBin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      angleBin = atof(tmpSs1);
      fAngleBinVector->push_back(angleBin);
      G4cout<<"angleBin = "<<angleBin/fAngleUnit<<G4endl;
    }
    else if ( position == tmpString.find("angleNo") ) // number of Tkin
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      angleNo = atoi(tmpSs1);
      jVector = 0;
      dataVector = new G4PhysicsFreeVector(angleNo);
      G4cout<<"angleNo = "<<angleNo<<G4endl;
    }
    else if ( position == tmpString.find("angle") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      angle = atof(tmpSs1)*fAngleUnit;
      fAngleVector->push_back(Tkin);
      G4cout<<"angle = "<<angle/degree<<G4endl;
    }
    */ ////////////////////////////////////////////////////////

    else if ( position == tmpString.find("no") ) // number of Tkin
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fNo = atoi(tmpSs1);
      kTable = 0;
      dataVector = new G4PhysicsFreeVector(fNo);
      G4cout<<"fNo = "<<fNo<<G4endl;
    }
    else if ( position == tmpString.find("Tkin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      Tkin = atof(tmpSs1);
      fTkinVector->push_back(Tkin);
      G4cout<<"Tkin = "<<Tkin<<G4endl;
    }
    else if ( tmpString.empty() ) // just skip empty lines
    {
      //     G4cout<<"tmpString is empty "<<G4endl;
      continue;
    }
    else // Tkin, ddxsc Filling to class field vectors
    { 
      G4String tmpStringone="";

      pos1 = tmpString.find('.');
      G4int postemp = tmpString.find(' ',pos1);
      tmpSs1   = tmpString.substr(0,postemp);
      tmpStringone=tmpString.substr(postemp,tmpString.length());

      //      G4cout << "Stringone" << tmpSs1 <<G4endl;

      pos2 = tmpStringone.find('.');
      postemp = tmpStringone.find(' ',pos2);
      tmpSs2   = tmpStringone.substr(0,postemp);
      tmpStringone=tmpStringone.substr(postemp,tmpStringone.length());
    
      //      G4cout << "Stringone" << tmpSs2 <<G4endl;
    
      pos3 = tmpStringone.find('.');
      postemp = tmpStringone.find(' ',pos3);
      tmpSs3   = tmpStringone.substr(0,postemp);
      tmpStringone=tmpStringone.substr(postemp,tmpStringone.length());

      //      G4cout << "Stringone" << tmpSs3 <<G4endl; 

      pos4 = tmpStringone.find('.');
      postemp = tmpStringone.find(' ',pos4);
      tmpSs4   = tmpStringone.substr(0,postemp);
      tmpStringone=tmpStringone.substr(postemp,tmpStringone.length());
      
      //      G4cout << "Stringone" << tmpSs4 <<G4endl; 
     
      tmpSs5   = tmpStringone;

      angle     = atof(tmpSs1);
      energy    = atof(tmpSs2);
      deltaE    = atof(tmpSs3);
      xsc       = atof(tmpSs4);
      deltaXsc  = atof(tmpSs5);

      // G4cout<<tmpSs1<<G4endl;

      //  G4cout<<tmpSs1<<"\t"<<tmpSs2<<"  "<<tmpSs3
      //     <<"  "<<tmpSs4<<"  "<<tmpSs5<<G4endl;

        G4cout<<angle<<"\t"<<energy<<"\t"<<deltaE
             <<"\t"<<xsc<<"\t"<<deltaXsc<<G4endl;

      fAngleVector->push_back(angle);
      fMomentumCVector->push_back(energy);
      fXscVector->push_back(xsc);       

      /*   
      xsc      = atof(tmpSs2);
        dataVector->PutValue(jVector,variable,xsc);
        jVector++;
        if(jVector == angleNo)
	{
          dataTable->push_back(dataVector);
          delete dataVector;
 
          kTable++;
          if( kTable == fNo )
	  {
            fDoubleDiffXscBank->push_back(dataTable);
            dataTable->clearAndDestroy();
	  }
	} 
	//  G4cout<<"variable = "<<variable<<"\t xsc = "<<xsc/millibarn<<G4endl;
      */
    }
  }
  fin.close();

  //sim Look if the file exists already and 
  // Append new data to already existing file

  struct stat s;
  int hfile;
  
  hfile=stat(fileName,&s);
  
    if(hfile !=-1)
      {
	G4cout<<"I found you"<< G4endl;
      }
  
  G4std::ofstream fout(fileName, G4std::ios::out | G4std::ios::app );

  // sim end

  fout.setf( G4std::ios::scientific, G4std::ios::floatfield );


  G4std::vector<size_t> deltaAngleVector;
  position = 0;

  angle = (*fAngleVector)[0];

  for(size_t i = 0; i < fAngleVector->size(); ++i)
  {
    if( angle != (*fAngleVector)[i] )
    {
      deltaAngleVector.push_back(i);
      angle = (*fAngleVector)[i];
    }
  }
  for( i = 0; i < fCommentVector->size(); ++i)
  {
    fout<<(*fCommentVector)[i]<<G4endl;
    G4cout<<(*fCommentVector)[i]<<G4endl;
  }

  fout<<G4endl<<"Tkin"<<"\t"<<(*fTkinVector)[0]<<G4endl;
  fout<<G4endl<<"no"<<"\t"<<deltaAngleVector.size()+1<<G4endl;



  fout<<G4endl<<"angleUnit"<<"\t"<<fAngleUnit<<G4endl;
  fout<<G4endl<<"energyUnit"<<"\t"<<fEnergyUnit<<G4endl;
  fout<<G4endl<<"ddXscUnit"<<"\t"<<fDdXscUnit<<G4endl;
  fout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[0]<<G4endl<<G4endl;
  fout<<G4endl<<"angleNo"<<"\t"<<deltaAngleVector[0]<<G4endl<<G4endl;

  G4cout<<G4endl<<"angleUnit"<<"\t"<<fAngleUnit<<G4endl;
  G4cout<<G4endl<<"energyUnit"<<"\t"<<fEnergyUnit<<G4endl;
  G4cout<<G4endl<<"ddXscUnit"<<"\t"<<fDdXscUnit<<G4endl;
  G4cout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[0]<<G4endl<<G4endl;
  G4cout<<G4endl<<"angleNo"<<"\t"<<deltaAngleVector[0]<<G4endl<<G4endl;

  i = 0;
  for(  size_t j = 0; j < fMomentumCVector->size(); ++j )
  {
    fout<<(*fMomentumCVector)[j]<<"\t"<<(*fXscVector)[j]<<G4endl;
    G4cout<<(*fMomentumCVector)[j]<<"\t"<<(*fXscVector)[j]<<G4endl;
    if( j == deltaAngleVector[i]-1 )
    { 
      i++;
      if (i < deltaAngleVector.size())
      {   
        fout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j+1]<<G4endl<<G4endl;
        G4cout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j+1]<<G4endl<<G4endl;
        fout<<G4endl<<"angleNo"<<"\t"<<deltaAngleVector[i]-deltaAngleVector[i-1]
           <<G4endl<<G4endl;
        G4cout<<G4endl<<"angleNo"<<"\t"<<deltaAngleVector[i]-deltaAngleVector[i-1]
          <<G4endl<<G4endl;
      }
      else
      {
        fout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j+1]<<G4endl<<G4endl;
        G4cout<<G4endl<<"angle"<<"\t"<<(*fAngleVector)[j+1]<<G4endl<<G4endl;
        fout<<G4endl<<"angleNo"<<"\t"<<fMomentumCVector->size()
                                          -deltaAngleVector[i-1]
           <<G4endl<<G4endl;
        G4cout<<G4endl<<"angleNo"<<"\t"<<fMomentumCVector->size()
                                           -deltaAngleVector[i-1]
          <<G4endl<<G4endl;
      }
    }
  }
  
  //sim Interpolation scheme 
  
  fout << G4endl << "# Space left for data interpolation scheme declaration"<< G4endl << G4endl;

  fout<<G4endl<<"# End data file" <<G4endl<<G4endl;

  /* 
  for(size_t i = 0; i < fDoubleDiffXscBank->size(); ++i)
  {
    for(  size_t j = 0; j < (*fDoubleDiffXscBank)[i]->length(); ++j )
    {
      for(size_t k = 0; k < (*(*fDoubleDiffXscBank)[i])(j)->GetVectorLength(); ++k)
      {
        G4cout<<(*fAngleVector)[i]/fAngleUnit<<"\t"
              <<(*(*fDoubleDiffXscBank)[i])(j)->GetLowEdgeEnergy(k)/fEnergyUnit
              <<"\t"
	      <<(*(*(*fDoubleDiffXscBank)[i])(j))(k)/fDdXscUnit<<G4endl;
      }
    } 
  }
  */
}

//////////////////////////////////////////////////////////

void G4VHadDataWriting::SetEnergyUnit(G4String& unitString)
{
  fEnergyUnit = unitString;
}

//////////////////////////////////////////////////////////

void G4VHadDataWriting::SetAngleUnit(G4String& unitString)
{
  fAngleUnit = unitString;
}

//////////////////////////////////////////////////////////

void G4VHadDataWriting::SetXscUnit(G4String& unitString)
{
  fXscUnit = unitString;
}

//////////////////////////////////////////////////////////

void G4VHadDataWriting::SetDdXscUnit(G4String& unitString)
{
  fDdXscUnit = unitString;
}



//
//
/////////////////////////////////////////////////////////////










