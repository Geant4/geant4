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
// $Id: HadrontherapyLet.cc
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), M. Sallemi, A. Salvia
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyLet.hh"
#include "HadrontherapyLetMessenger.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"

#define D 1.18 //PMMA density
#define N 200 //number of slices in which phantom is divided

G4double vetfluence[N];

/////////////////////////////////////////////////////////////////////////////
HadrontherapyLet::HadrontherapyLet(HadrontherapyPrimaryGeneratorAction* PGApointer)
{
  pga=PGApointer;  //pointer to H.PrimaryGeneratorAction class
  letMessenger = new HadrontherapyLetMessenger(this);
  
  ffluence.open("fluence.txt"); 

  for(i=0;i<N;i++)   vetfluence[i]=0;
 
}
/////////////////////////////////////////////////////////////////////////////
HadrontherapyLet::~HadrontherapyLet()
{
  ffluence.close();

  for (i=0;i<size;i++)  fspectrum[i].close();
  delete fspectrum;
     
  for(i=0;i<size;i++) delete spectrum[i];
  delete spectrum;

}

/////////////////////////////////////////////////////////////////////////////
void  HadrontherapyLet::Fluence_Let(G4int x, G4double energy) {

  vetfluence[x]++;  //fluence calculation

  for(i=0;i<size;i++)
   {
    if ( x==((int)(vetdepth.at(i)/0.2)) )    
         { fspectrum[i] << energy << G4endl;
	   j=(int)(energy/0.25);
	   spectrum[i][j]++;
	   return;
         }
   }
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyLet::LetOutput()
{
  for(i=0;i<N;i++)
    ffluence << i << '\t' << vetfluence[i] << G4endl; 

  G4cout << "File di output fluence.txt creato."<<G4endl;
  
 if(size)
 {
  
  fstop.open("stopPMMA.in"); //"stopWATER.in" for stopping power in water (change density!)
                             //stopping power file corresponding to mean energy of 0.25 MeV bins

  stop = new G4double[bins];

  for(i=0;i<bins;i++)
    {
        fstop>>stop[i];      //read stopping power from file
	stop[i]=stop[i]*D;   //multiply stopping power for density                    
    }
  fstop.close();

  //LET calculation

  for(i=0;i<size;i++)
    {
      n1=0;
      d1=0;
      n2=0;
      d2=0;
      for(j=0;j<bins;j++)
      {
        n1=n1+(spectrum[i][j]*stop[j]);
        d1=d1+spectrum[i][j];
	
	n2=n2+(spectrum[i][j]*stop[j]*stop[j]);
	d2=d2+(spectrum[i][j]*stop[j]);
	  
       } 
 
      lett=(n1/d1)/10;  // divide let per 10 to obtain KeV/micrometers
      letd=(n2/d2)/10;  // instead of MeV/cm (obtained from stopping power)
     flet << lett <<'\t' << letd << G4endl;
     
    }
     flet.close(); 
     G4cout << "File di output let.txt creato." << G4endl;
 }

}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyLet::SetValue(G4String value)
{
  G4String str(value),nums;
  str_size k=0,l;
  G4double num;

  ostringstream ost;

  flet.open("let.txt");  //controllo su apertura o esistenza file!
  
  while (str.index("-",k)!=std::string::npos)
  {
  l=str.index("-",k);

  nums=str(k,l-k);
  k=l+1;
  
  stringstream sst; 
  sst<<nums;
  sst>>num;
  vetdepth.push_back(num);
  }

  nums=str(k,str.length());
 
 stringstream sst;
 sst<<nums;
 sst>>num; 
 vetdepth.push_back(num);

 size=(int)vetdepth.size();
 
 G4double defenergy = pga -> GetmeanKineticEnergy();
 G4int energylimit =((G4int)(defenergy/10)+1)*10;

 bins = (G4int)(energylimit/0.25);
  
//allocazione dinamica vettore di spettri
 spectrum = new G4int*[size];

 for(i=0;i<size;i++)
  {
   spectrum[i]=new G4int[bins];
   for(j=0;j<bins;j++)
     spectrum[i][j]=0;
  }

 //allocazione dinamica vettore di file

 fspectrum = new ofstream[size];
  
 for(i=1;i<=size;i++)
    { ost << i;
      nome_file="spettro"+ ost.str() + ".txt";
      ost.str("");
      fspectrum[i-1].open(nome_file);
    }

}

