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

#include "g4std/fstream"
#include "g4std/strstream"
#include "globals.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


G4HadDataReading::G4HadDataReading()
{
  fEnergyUnit     = GeV; 
  fAngleUnit      = degree; 
  fXscUnit        = millibarn;
  fXscPeAngleUnit = millibarn/degree; 
  fXscPerMomCUnit = millibarn/GeV;
  fDdXscUnit      = millibarn/degree/GeV;

  fTkinBin = -1.;
  fNo = -1;
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
}



G4HadDataReading::~G4HadDataReading()
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
  G4cout<<"G4HadDataReading::LoadDifferentialXSC(fN),\n"<<" fN = "
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
  G4double Tkin, variable, xsc, angleBin, angle;
  G4int kTable, jVector, angleNo;
  G4PhysicsFreeVector* dataVector = NULL;
  G4PhysicsTable*  dataTable = NULL;

  while( !fin.eof() )
  {
    getline(fin,tmpString);
    G4cout<<tmpString<<G4endl;
    
    if( position == tmpString.find('#') ||
        position == tmpString.find('*') ) 
    {
      G4cout<<tmpString<<G4endl;  // print comment lines
    }
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
    else if ( position == tmpString.find("no") ) // number of Tkin
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      fNo = atoi(tmpSs1);
      kTable = 0;
      dataTable = new G4PhysicsTable(fNo);
      G4cout<<"fNo = "<<fNo<<G4endl;
    }
    else if ( position == tmpString.find("Tkin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      Tkin = atof(tmpSs1)*fEnergyUnit;
      fTkinVector->push_back(Tkin);
      G4cout<<"Tkin = "<<Tkin/fEnergyUnit<<G4endl;
    }
    // Treatment of different angles at a given energy

    else if ( position == tmpString.find("angleBin") )
    {
      posValue = tmpString.find("\t");
      tmpSs1   = tmpString.substr(posValue+1);
      angleBin = atof(tmpSs1)*fAngleUnit;
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
      fAngleVector->push_back(angle);
      G4cout<<"angle = "<<angle/degree<<G4endl;
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

        variable     = atof(tmpSs1)*fEnergyUnit;
        xsc      = atof(tmpSs2)*fDdXscUnit;
        dataVector->PutValue(jVector,variable,xsc);
        jVector++;
        if(jVector == angleNo)
	{
          dataTable->push_back(dataVector);
	  //  delete dataVector;
          jVector = 0;
          kTable++;

	  //  if( kTable == fNo )
          if( (G4int)dataTable->size() == fNo )
	  {
            fDoubleDiffXscBank->push_back(dataTable);
	    //  dataTable->clearAndDestroy();
            kTable = 0;
	  }
	} 
	//  G4cout<<"variable = "<<variable<<"\t xsc = "<<xsc/millibarn<<G4endl;
        
    }
  }
  fin.close();

  G4cout<<"Output of fDoubleDiffXscBank data"<<G4endl;
  G4cout<<"fDoubleDiffXscBank->size() = "<<fDoubleDiffXscBank->size()<<G4endl;
 
  for(size_t i = 0; i < fDoubleDiffXscBank->size(); ++i)
  {
    G4cout<<"(*fDoubleDiffXscBank)["<<i<<"]->size() = "
        <<(*fDoubleDiffXscBank)[i]->size()<<G4endl;

    for(  size_t j = 0; j < (*fDoubleDiffXscBank)[i]->size(); ++j )
    {
      G4cout<<"(*(*fDoubleDiffXscBank)["<<i<<"])("<<j<<")->GetVectorLength() = "
        <<(*(*fDoubleDiffXscBank)[i])(j)->GetVectorLength()<<G4endl;

      for(size_t k = 0; k < (*(*fDoubleDiffXscBank)[i])(j)->GetVectorLength(); ++k)
      {
        G4cout<<(*fAngleVector)[j]/fAngleUnit<<"\t"
              <<(*(*fDoubleDiffXscBank)[i])(j)->GetLowEdgeEnergy(k)/fEnergyUnit
              <<"\t"
	      <<(*(*(*fDoubleDiffXscBank)[i])(j))(k)/fDdXscUnit<<G4endl;
      }
    } 
  }
}

//////////////////////////////////////////////////////////

void G4HadDataReading::SetEnergyUnit(G4String& unitString)
{
  if      ( unitString == "MeV" ) fEnergyUnit = MeV;
  else if ( unitString == "keV" ) fEnergyUnit = keV;
  else if ( unitString == "eV" )  fEnergyUnit = eV;
  else if ( unitString == "TeV" ) fEnergyUnit = TeV;
  else if ( unitString == "PeV" ) fEnergyUnit = PeV;
}

//////////////////////////////////////////////////////////

void G4HadDataReading::SetAngleUnit(G4String& unitString)
{
  if      ( unitString == "rad" )  fAngleUnit = rad;
  else if ( unitString == "mrad" ) fAngleUnit = mrad;
  else if ( unitString == "sr" )   fAngleUnit = sr;
  else if ( unitString == "deg" )  fAngleUnit = degree;
}

//////////////////////////////////////////////////////////

void G4HadDataReading::SetXscUnit(G4String& unitString)
{
  if      ( unitString == "barn" )      fXscUnit = barn;
  else if ( unitString == "microbarn" ) fXscUnit = microbarn;
  else if ( unitString == "nanobarn" )  fXscUnit = nanobarn;
  else if ( unitString == "picobarn" )  fXscUnit = picobarn;
}

//////////////////////////////////////////////////////////

void G4HadDataReading::SetDdXscUnit(G4String& unitString)
{
  if      ( unitString == "barn/MeV/sr" )      fXscUnit = barn/MeV/sr;
  else if ( unitString == "microbarn/MeV/sr" ) fXscUnit = microbarn/MeV/sr;
  else if ( unitString == "nanobarn/MeV/sr" )  fXscUnit = nanobarn/MeV/sr;
  else if ( unitString == "picobarn/MeV/sr" )  fXscUnit = picobarn/MeV/sr;
}



//
//
/////////////////////////////////////////////////////////////










