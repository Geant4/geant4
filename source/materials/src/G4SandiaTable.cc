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

// $Id: G4SandiaTable.cc,v 1.43 2010-12-23 16:12:55 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// 10.06.97 created. V. Grichine
// 18.11.98 simplified public interface; new methods for materials.  mma
// 31.01.01 redesign of ComputeMatSandiaMatrix().  mma
// 16.02.01 adapted for STL.  mma
// 22.02.01 GetsandiaCofForMaterial(energy) return 0 below lowest interval  mma  
// 03.04.01 fnulcof returned if energy < emin
// 10.07.01 Migration to STL. M. Verderi.
// 03.02.04 Update distructor V.Ivanchenko
// 05.03.04 New methods for old sorting algorithm for PAI model. V.Grichine
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....


#include "G4SandiaTable.hh"
#include "G4StaticSandiaData.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

G4int    G4SandiaTable::fCumulInterval[101]  = {0};
G4double G4SandiaTable::fSandiaCofPerAtom[4] = {0.0};
G4double const G4SandiaTable::funitc[4] = {cm2*keV/g,     
					   cm2*keV*keV/g,     
					   cm2*keV*keV*keV/g,     
					   cm2*keV*keV*keV*keV/g};
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4SandiaTable::G4SandiaTable(G4Material* material)
  : fMaterial(material)
{
  fMatSandiaMatrix    = 0; 
  fMatSandiaMatrixPAI = 0;
  fPhotoAbsorptionCof = 0;

  fMatNbOfIntervals   = 0;

 
  fMaxInterval        = 0;
  fVerbose            = 0;  

  //build the CumulInterval array

  fCumulInterval[0] = 1;

  for (G4int Z=1; Z<101; ++Z) {
    fCumulInterval[Z] = fCumulInterval[Z-1] + fNbOfIntervals[Z];
  }
  
  //initialisation of fnulcof
  fnulcof[0] = fnulcof[1] = fnulcof[2] = fnulcof[3] = 0.;

  fMaxInterval = 0;

  //compute macroscopic Sandia coefs for a material   
  ComputeMatSandiaMatrix(); // mma

  // ComputeMatTable();  // vmg
}
							
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4SandiaTable::G4SandiaTable(__void__&)
  : fMaterial(0),fMatSandiaMatrix(0),fMatSandiaMatrixPAI(0),fPhotoAbsorptionCof(0)
{
  fnulcof[0] = fnulcof[1] = fnulcof[2] = fnulcof[3] = 0.;
  fMaxInterval = 0;
  fMatNbOfIntervals = 0;
  fVerbose          = 0;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4SandiaTable::~G4SandiaTable()
{ 
  if(fMatSandiaMatrix) {
    fMatSandiaMatrix->clearAndDestroy();
    delete fMatSandiaMatrix;
  }
  if(fMatSandiaMatrixPAI) {
    fMatSandiaMatrixPAI->clearAndDestroy();
    delete fMatSandiaMatrixPAI;
  }
  if(fPhotoAbsorptionCof)
  {
    delete [] fPhotoAbsorptionCof;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double G4SandiaTable::GetZtoA(G4int Z)
{
  return fZtoAratio[Z];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double*
G4SandiaTable::GetSandiaCofPerAtom(G4int Z, G4double energy)
{
  assert (Z > 0 && Z < 101);
   
  G4double Emin  = fSandiaTable[fCumulInterval[Z-1]][0]*keV;
  G4double Iopot = fIonizationPotentials[Z]*eV;
  if (Iopot > Emin) Emin = Iopot;
   
  G4int interval = fNbOfIntervals[Z] - 1;
  G4int row = fCumulInterval[Z-1] + interval;
  while ((interval>0) && (energy<fSandiaTable[row][0]*keV)) {
    --interval;
    row = fCumulInterval[Z-1] + interval;
  }
  if (energy >= Emin)
    {        
      G4double AoverAvo = Z*amu/fZtoAratio[Z];
         
      fSandiaCofPerAtom[0]=AoverAvo*funitc[0]*fSandiaTable[row][1];     
      fSandiaCofPerAtom[1]=AoverAvo*funitc[1]*fSandiaTable[row][2];     
      fSandiaCofPerAtom[2]=AoverAvo*funitc[2]*fSandiaTable[row][3];     
      fSandiaCofPerAtom[3]=AoverAvo*funitc[3]*fSandiaTable[row][4];
    }
  else 
    {
      fSandiaCofPerAtom[0] = fSandiaCofPerAtom[1] = fSandiaCofPerAtom[2] =
	fSandiaCofPerAtom[3] = 0.;
    }                
  return fSandiaCofPerAtom;     
}
						 	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4SandiaTable::ComputeMatSandiaMatrix()
{  
  //get list of elements
  const G4int NbElm = fMaterial->GetNumberOfElements();
  const G4ElementVector* ElementVector = fMaterial->GetElementVector();
  
  G4int* Z = new G4int[NbElm];               //Atomic number
   
  //determine the maximum number of energy-intervals for this material
  
  G4int MaxIntervals = 0;
  G4int elm;
    
  for ( elm = 0; elm < NbElm; elm++ ) 
  {
    Z[elm] = (G4int)(*ElementVector)[elm]->GetZ();
    MaxIntervals += fNbOfIntervals[Z[elm]];
  }  
     
  //copy the Energy bins in a tmp1 array
  //(take care of the Ionization Potential of each element)
  
  G4double* tmp1 = new G4double[MaxIntervals]; 
  G4double IonizationPot;
  G4int interval1 = 0;

  for ( elm = 0; elm < NbElm; elm++ ) 
  {
    IonizationPot = GetIonizationPot(Z[elm]);

    for (G4int row = fCumulInterval[ Z[elm]-1 ]; row < fCumulInterval[Z[elm]]; row++ ) 
    {
      tmp1[interval1++] = std::max(fSandiaTable[row][0]*keV,IonizationPot);
    }
  }   
        
  //sort the energies in strickly increasing values in a tmp2 array
  //(eliminate redondances)
  
  G4double* tmp2 = new G4double[MaxIntervals];
  G4double Emin;
  G4int interval2 = 0;
  
  do 
  {
    Emin = DBL_MAX;

    for ( G4int i1 = 0; i1 < MaxIntervals; i1++ ) 
    {
      if (tmp1[i1] < Emin) Emin = tmp1[i1];          //find the minimum
    }
    if (Emin < DBL_MAX) tmp2[interval2++] = Emin;
    //copy Emin in tmp2
    for ( G4int j1 = 0; j1 < MaxIntervals; j1++ ) 
    {
      if (tmp1[j1] <= Emin) tmp1[j1] = DBL_MAX;      //eliminate from tmp1	    
    }
  } while (Emin < DBL_MAX);
        	         	   	
  //create the sandia matrix for this material
    
  fMatSandiaMatrix = new G4OrderedTable();
  G4int interval;

  for (interval = 0; interval < interval2; interval++ ) 
  {
    fMatSandiaMatrix->push_back( new G4DataVector(5,0.) );
  }
  
  //ready to compute the Sandia coefs for the material
  
  const G4double* NbOfAtomsPerVolume = fMaterial->GetVecNbOfAtomsPerVolume();
  
  const G4double prec = 1.e-03*eV;
  G4double coef, oldsum(0.), newsum(0.);
  fMatNbOfIntervals = 0;
         
  for ( interval = 0; interval < interval2; interval++ ) 
  {
    Emin = (*(*fMatSandiaMatrix)[fMatNbOfIntervals])[0] = tmp2[interval];

    for ( G4int k = 1; k < 5; k++ )   (*(*fMatSandiaMatrix)[fMatNbOfIntervals])[k] = 0.;      
    
    newsum = 0.;
      
    for ( elm = 0; elm < NbElm; elm++ ) 
    {    
      GetSandiaCofPerAtom(Z[elm], Emin+prec);

      for ( G4int j = 1; j < 5; j++ ) 
      {
	coef = NbOfAtomsPerVolume[elm]*fSandiaCofPerAtom[j-1];
	(*(*fMatSandiaMatrix)[fMatNbOfIntervals])[j] += coef;
	newsum += std::abs(coef);
      }						       
    }	      			      			      	 
    //check for null or redondant intervals
	 
    if (newsum != oldsum) { oldsum = newsum; fMatNbOfIntervals++;}
  }
  delete [] Z;
  delete [] tmp1;
  delete [] tmp2;

  if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
  {
    G4cout<<"mma, G4SandiaTable::ComputeMatSandiaMatrix(), mat = "<<fMaterial->GetName()<<G4endl;

    for( G4int i = 0; i < fMatNbOfIntervals; i++)
    {
      G4cout<<i<<"\t"<<GetSandiaCofForMaterial(i,0)/keV<<" keV \t"<<this->GetSandiaCofForMaterial(i,1)
       <<"\t"<<this->GetSandiaCofForMaterial(i,2)<<"\t"<<this->GetSandiaCofForMaterial(i,3)
       <<"\t"<<this->GetSandiaCofForMaterial(i,4)<<G4endl;
    }   
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4SandiaTable::ComputeMatSandiaMatrixPAI()
{  
  G4int MaxIntervals = 0;
  G4int elm, c, i, j, jj, k, k1, k2, c1, n1;    

  const G4int noElm = fMaterial->GetNumberOfElements();
  const G4ElementVector* ElementVector = fMaterial->GetElementVector();  
  G4int* Z = new G4int[noElm];               //Atomic number

  for ( elm = 0; elm < noElm; elm++ )
  { 
    Z[elm] = (G4int)(*ElementVector)[elm]->GetZ();
    MaxIntervals += fNbOfIntervals[Z[elm]];
  }  
  fMaxInterval = MaxIntervals + 2;

  if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
  {
    G4cout<<"fMaxInterval = "<<fMaxInterval<<G4endl;
  } 
  G4double* fPhotoAbsorptionCof0 = new G4double[fMaxInterval];
  G4double* fPhotoAbsorptionCof1 = new G4double[fMaxInterval];
  G4double* fPhotoAbsorptionCof2 = new G4double[fMaxInterval];
  G4double* fPhotoAbsorptionCof3 = new G4double[fMaxInterval];
  G4double* fPhotoAbsorptionCof4 = new G4double[fMaxInterval];

  for( c = 0; c < fMaxInterval; c++ )   // just in case
  {
    fPhotoAbsorptionCof0[c] = 0.;
    fPhotoAbsorptionCof1[c] = 0.;
    fPhotoAbsorptionCof2[c] = 0.;
    fPhotoAbsorptionCof3[c] = 0.;
    fPhotoAbsorptionCof4[c] = 0.;
  }
  c = 1;

  for(i = 0; i < noElm; i++)
  {
    G4double I1 = fIonizationPotentials[Z[i]]*keV;  // I1 in keV
    n1          = 1;                                    
 
    for( j = 1; j < Z[i]; j++ )  n1 += fNbOfIntervals[j];
    
    G4int n2 = n1 + fNbOfIntervals[Z[i]];
    
    for( k1 = n1; k1 < n2; k1++ )
    {
      if( I1  > fSandiaTable[k1][0] )
      {
	 continue;    // no ionization for energies smaller than I1 (first
      }		       // ionisation potential)		     
      break;
    }
    G4int flag = 0;
    
    for( c1 = 1; c1 < c; c1++ )
    {
      if( fPhotoAbsorptionCof0[c1] == I1 ) // this value already has existed
      {
	flag = 1;                      
	break;                         
      }
    }
    if(flag == 0)
    {
      fPhotoAbsorptionCof0[c] = I1;
      c++;
    }
    for( k2 = k1; k2 < n2; k2++ )
    {
      flag = 0;

      for( c1 = 1; c1 < c; c1++ )
      {
        if( fPhotoAbsorptionCof0[c1] == fSandiaTable[k2][0] )
        {
	  flag = 1;
 	  break;
        }
      }
      if(flag == 0)
      {
        fPhotoAbsorptionCof0[c] = fSandiaTable[k2][0];
	c++;
      }
    }       
  }   // end for(i)
  // sort out

  for( i = 1; i < c; i++ ) 
  {
    for( j = i + 1; j < c;  j++ )
    {
      if( fPhotoAbsorptionCof0[i] > fPhotoAbsorptionCof0[j] ) 
      {
        G4double tmp = fPhotoAbsorptionCof0[i];
	fPhotoAbsorptionCof0[i] = fPhotoAbsorptionCof0[j];
	fPhotoAbsorptionCof0[j] = tmp;
      }
    }
    if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
    {
      G4cout<<i<<"\t energy = "<<fPhotoAbsorptionCof0[i]<<G4endl;
    }
  } 
  fMaxInterval = c;
 
  const G4double* fractionW = fMaterial->GetFractionVector();

  if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
  {
    for( i = 0; i < noElm; i++ )  G4cout<<i<<" = elN, fraction = "<<fractionW[i]<<G4endl;
  }      
   
  for( i = 0; i < noElm; i++ )
  {
    n1 = 1;
    G4double I1 = fIonizationPotentials[Z[i]]*keV;

    for( j = 1; j < Z[i]; j++ )  n1 += fNbOfIntervals[j];
	
    G4int n2 = n1 + fNbOfIntervals[Z[i]] - 1;

    for(k = n1; k < n2; k++)
    {
      G4double B1 = fSandiaTable[k][0];
      G4double B2 = fSandiaTable[k+1][0];

      for(G4int c = 1; c < fMaxInterval-1; c++)
      {
	G4double E1 = fPhotoAbsorptionCof0[c];
	G4double E2 = fPhotoAbsorptionCof0[c+1];

        if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
        {
          G4cout<<"k = "<<k<<", c = "<<c<<", B1 = "<<B1<<", B2 = "<<B2<<", E1 = "<<E1<<", E2 = "<<E2<<G4endl;
        }      
	if( B1 > E1 || B2 < E2 || E1 < I1 )  
	{
          if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
          {
            G4cout<<"continue for: B1 = "<<B1<<", B2 = "<<B2<<", E1 = "<<E1<<", E2 = "<<E2<<G4endl;
          }      
          continue;
	}		
	fPhotoAbsorptionCof1[c] += fSandiaTable[k][1]*fractionW[i];
	fPhotoAbsorptionCof2[c] += fSandiaTable[k][2]*fractionW[i];
	fPhotoAbsorptionCof3[c] += fSandiaTable[k][3]*fractionW[i];
	fPhotoAbsorptionCof4[c] += fSandiaTable[k][4]*fractionW[i];
      }  
    }   
    // Last interval

    fPhotoAbsorptionCof1[fMaxInterval-1] += fSandiaTable[k][1]*fractionW[i];
    fPhotoAbsorptionCof2[fMaxInterval-1] += fSandiaTable[k][2]*fractionW[i];
    fPhotoAbsorptionCof3[fMaxInterval-1] += fSandiaTable[k][3]*fractionW[i];
    fPhotoAbsorptionCof4[fMaxInterval-1] += fSandiaTable[k][4]*fractionW[i];
  }     // for(i)  
  c = 0;     // Deleting of first intervals where all coefficients = 0

  do                        
  {
    c++;

    if( fPhotoAbsorptionCof1[c] != 0.0 ||
	fPhotoAbsorptionCof2[c] != 0.0 ||
        fPhotoAbsorptionCof3[c] != 0.0 || 
	fPhotoAbsorptionCof4[c] != 0.0     )  continue;

    if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
    {
      G4cout<<c<<" = number with zero cofs"<<G4endl;
    }      
    for( jj = 2; jj < fMaxInterval; jj++ )
    {
      fPhotoAbsorptionCof0[jj-1] = fPhotoAbsorptionCof0[jj];
      fPhotoAbsorptionCof1[jj-1] = fPhotoAbsorptionCof1[jj];
      fPhotoAbsorptionCof2[jj-1] = fPhotoAbsorptionCof2[jj];
      fPhotoAbsorptionCof3[jj-1] = fPhotoAbsorptionCof3[jj];
      fPhotoAbsorptionCof4[jj-1] = fPhotoAbsorptionCof4[jj];
    }
    fMaxInterval--;
    c--;
  }
  while( c < fMaxInterval - 1 );
  	
  // create the sandia matrix for this material
    
  fMatSandiaMatrixPAI = new G4OrderedTable();
  G4double density = fMaterial->GetDensity();
 
  for (i = 0; i < fMaxInterval; i++)  fMatSandiaMatrixPAI->push_back(new G4DataVector(5,0.));
    	         	
  for (i = 0; i < fMaxInterval; i++)
  {
    (*(*fMatSandiaMatrixPAI)[i])[0] = fPhotoAbsorptionCof0[i+1];
    (*(*fMatSandiaMatrixPAI)[i])[1] = fPhotoAbsorptionCof1[i+1]*density;
    (*(*fMatSandiaMatrixPAI)[i])[2] = fPhotoAbsorptionCof2[i+1]*density;
    (*(*fMatSandiaMatrixPAI)[i])[3] = fPhotoAbsorptionCof3[i+1]*density;
    (*(*fMatSandiaMatrixPAI)[i])[4] = fPhotoAbsorptionCof4[i+1]*density;
  }
  if ( fVerbose > 0 && fMaterial->GetName() == "G4_Ar" )
  {
    G4cout<<"mma, G4SandiaTable::ComputeMatSandiaMatrixPAI(), mat = "<<fMaterial->GetName()<<G4endl;

    for( G4int i = 0; i < fMaxInterval; i++)
    {
      G4cout<<i<<"\t"<<GetSandiaMatTablePAI(i,0)/keV<<" keV \t"<<this->GetSandiaMatTablePAI(i,1)
       <<"\t"<<this->GetSandiaMatTablePAI(i,2)<<"\t"<<this->GetSandiaMatTablePAI(i,3)
       <<"\t"<<this->GetSandiaMatTablePAI(i,4)<<G4endl;
    }   
  }

	         	    
  delete [] Z;
  delete [] fPhotoAbsorptionCof0;
  delete [] fPhotoAbsorptionCof1;
  delete [] fPhotoAbsorptionCof2;
  delete [] fPhotoAbsorptionCof3;
  delete [] fPhotoAbsorptionCof4;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//
// Methods for PAI model

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4SandiaTable::G4SandiaTable(G4int matIndex)
{ 
  fMaterial           = 0;
  fMatNbOfIntervals   = 0;
  fMatSandiaMatrix    = 0; 
  fMatSandiaMatrixPAI = 0;
  fPhotoAbsorptionCof = 0;

 
  fMaxInterval        = 0;
  fVerbose            = 0;  

  //initialisation of fnulcof
  fnulcof[0] = fnulcof[1] = fnulcof[2] = fnulcof[3] = 0.;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int numberOfMat = G4Material::GetNumberOfMaterials();

  if ( matIndex >= 0 && matIndex < numberOfMat)
    {
      fMaterial = (*theMaterialTable)[matIndex];
      ComputeMatTable();
    }
  else
    {
      G4Exception("G4SandiaTable::G4SandiaTable(G4int matIndex): wrong matIndex ");
    }
}

///////////////////////////////////////////////////////////////////////
//
// Bubble sorting of left energy interval in SandiaTable in ascening order
//

void
G4SandiaTable::SandiaSort(G4double** da ,
 			  G4int sz )
{
  for(G4int i = 1;i < sz; i++ ) 
   {
     for(G4int j = i + 1;j < sz; j++ )
     {
       if(da[i][0] > da[j][0])   SandiaSwap(da,i,j);      
     }
   }
}

////////////////////////////////////////////////////////////////////////////
//
//  SandiaIntervals 
//

G4int
G4SandiaTable::SandiaIntervals(G4int Z[],
			       G4int el )
{
  G4int c,  i, flag = 0, n1 = 1;
  G4int j, c1, k1, k2;
  G4double I1;
  fMaxInterval = 0;

  for( i = 0; i < el; i++ )  fMaxInterval += fNbOfIntervals[ Z[i] ]; 

  fMaxInterval += 2;

  if( fVerbose > 0 ) G4cout<<"begin sanInt, fMaxInterval = "<<fMaxInterval<<G4endl;

  fPhotoAbsorptionCof = new G4double* [fMaxInterval];

  for( i = 0; i < fMaxInterval; i++ )   fPhotoAbsorptionCof[i] = new G4double[5];
 
  //  for(c = 0; c < fIntervalLimit; c++)   // just in case

  for( c = 0; c < fMaxInterval; c++ )     fPhotoAbsorptionCof[c][0] = 0.;
  
  c = 1;

  for( i = 0; i < el; i++ )
  {
    I1 = fIonizationPotentials[ Z[i] ]*keV;  // First ionization
    n1 = 1;                                  // potential in keV

    for( j = 1; j < Z[i]; j++ )  n1 += fNbOfIntervals[j];
    
    G4int n2 = n1 + fNbOfIntervals[Z[i]];
    
    for( k1 = n1; k1 < n2; k1++ )
    {
      if( I1  > fSandiaTable[k1][0] )
      {
	 continue;    // no ionization for energies smaller than I1 (first
      }		       // ionisation potential)		     
      break;
    }
    flag = 0;
    
    for( c1 = 1; c1 < c; c1++ )
    {
      if( fPhotoAbsorptionCof[c1][0] == I1 ) // this value already has existed
      {
	flag = 1;                      
	break;                         
      }
    }
    if( flag == 0 )
    {
      fPhotoAbsorptionCof[c][0] = I1;
      c++;
    }
    for( k2 = k1; k2 < n2; k2++ )
    {
      flag = 0;

      for( c1 = 1; c1 < c; c1++ )
      {
        if( fPhotoAbsorptionCof[c1][0] == fSandiaTable[k2][0] )
        {
	  flag = 1;
 	  break;
        }
      }
      if( flag == 0 )
      {
        fPhotoAbsorptionCof[c][0] = fSandiaTable[k2][0];
	if( fVerbose > 0 ) G4cout<<"sanInt, c = "<<c<<", E_c = "<<fPhotoAbsorptionCof[c][0]<<G4endl;
	c++;
      }
    }       
  }   // end for(i)
  
  SandiaSort(fPhotoAbsorptionCof,c);
  fMaxInterval = c;
  if( fVerbose > 0 ) G4cout<<"end SanInt, fMaxInterval = "<<fMaxInterval<<G4endl;
  return c;
}   

///////////////////////////////////////////////////////////////////////
//
//  SandiaMixing
//

G4int
G4SandiaTable::SandiaMixing(         G4int Z[],
			       const G4double fractionW[],
			             G4int el,
			             G4int mi     )
{
  G4int i, j, n1, k, c=1, jj, kk;
  G4double I1, B1, B2, E1, E2;
   
  for( i = 0; i < mi; i++ )
  {
    for( j = 1; j < 5; j++ ) fPhotoAbsorptionCof[i][j] = 0.;
  }
  for( i = 0; i < el; i++ )
  {
    n1 = 1;
    I1 = fIonizationPotentials[Z[i]]*keV;

    for( j = 1; j < Z[i]; j++ )   n1 += fNbOfIntervals[j];
      
    G4int n2 = n1 + fNbOfIntervals[Z[i]] - 1;

    for( k = n1; k < n2; k++ )
    {
      B1 = fSandiaTable[k][0];
      B2 = fSandiaTable[k+1][0];

      for( c = 1; c < mi-1; c++ )
      {
        E1 = fPhotoAbsorptionCof[c][0];
        E2 = fPhotoAbsorptionCof[c+1][0];

        if( B1 > E1 || B2 < E2 || E1 < I1 )   continue;
	    
	for( j = 1; j < 5; j++ ) 
	{
          fPhotoAbsorptionCof[c][j] += fSandiaTable[k][j]*fractionW[i];
          if( fVerbose > 0 )
	  {
	    G4cout<<"c="<<c<<"; j="<<j<<"; fST="<<fSandiaTable[k][j]<<"; frW="<<fractionW[i]<<G4endl;
	  }
	}	    
      }  
    }   
    for( j = 1; j < 5; j++ )   // Last interval
    {
      fPhotoAbsorptionCof[mi-1][j] += fSandiaTable[k][j]*fractionW[i];
      if( fVerbose > 0 )
      {
	G4cout<<"mi-1="<<mi-1<<"; j="<<j<<"; fST="<<fSandiaTable[k][j]<<"; frW="<<fractionW[i]<<G4endl;
      }
    }
  }     // for(i)
  c = 0;     // Deleting of first intervals where all coefficients = 0

  do                        
  {
    c++;

    if( fPhotoAbsorptionCof[c][1] != 0.0 ||
	fPhotoAbsorptionCof[c][2] != 0.0 ||
        fPhotoAbsorptionCof[c][3] != 0.0 || 
	fPhotoAbsorptionCof[c][4] != 0.0     )  continue;
       
    for( jj = 2; jj < mi; jj++ )
    {
      for( kk = 0; kk < 5; kk++ ) fPhotoAbsorptionCof[jj-1][kk] = fPhotoAbsorptionCof[jj][kk];	  
    }
    mi--;
    c--;
  }
  while( c < mi - 1 );

  if( fVerbose > 0 ) G4cout<<"end SanMix, mi = "<<mi<<G4endl;
    
  return mi;
}  

////////////////////////////////////////////////////////////////////////////
//
//  Sandia interval and mixing calculations for materialCutsCouple constructor 
//

void G4SandiaTable::ComputeMatTable()
{
  G4int MaxIntervals = 0;
  G4int elm, c, i, j, jj, k, kk, k1, k2, c1, n1;    

  const G4int noElm = fMaterial->GetNumberOfElements();
  const G4ElementVector* ElementVector = fMaterial->GetElementVector();  
  G4int* Z = new G4int[noElm];               //Atomic number

  for (elm = 0; elm<noElm; elm++)
  { 
    Z[elm] = (G4int)(*ElementVector)[elm]->GetZ();
    MaxIntervals += fNbOfIntervals[Z[elm]];
  }  
  fMaxInterval = 0;

  for(i = 0; i < noElm; i++)  fMaxInterval += fNbOfIntervals[Z[i]]; 
  
  fMaxInterval += 2;

//  G4cout<<"fMaxInterval = "<<fMaxInterval<<G4endl;

  fPhotoAbsorptionCof = new G4double* [fMaxInterval];

  for(i = 0; i < fMaxInterval; i++)  
  {
     fPhotoAbsorptionCof[i] = new G4double[5];
  }

  //  for(c = 0; c < fIntervalLimit; c++)   // just in case

  for(c = 0; c < fMaxInterval; c++)   // just in case
  {
     fPhotoAbsorptionCof[c][0] = 0.;
  }
  c = 1;

  for(i = 0; i < noElm; i++)
  {
    G4double I1 = fIonizationPotentials[Z[i]]*keV;  // First ionization
    n1 = 1;                                     // potential in keV
 
    for(j = 1; j < Z[i]; j++)
    {
      n1 += fNbOfIntervals[j];
    }
    G4int n2 = n1 + fNbOfIntervals[Z[i]];
    
    for(k1 = n1; k1 < n2; k1++)
    {
      if(I1  > fSandiaTable[k1][0])
      {
	 continue;    // no ionization for energies smaller than I1 (first
      }		       // ionisation potential)		     
      break;
    }
    G4int flag = 0;
    
    for(c1 = 1; c1 < c; c1++)
    {
      if(fPhotoAbsorptionCof[c1][0] == I1) // this value already has existed
      {
	flag = 1;                      
	break;                         
      }
    }
    if(flag == 0)
    {
      fPhotoAbsorptionCof[c][0] = I1;
      c++;
    }
    for(k2 = k1; k2 < n2; k2++)
    {
      flag = 0;

      for(c1 = 1; c1 < c; c1++)
      {
        if(fPhotoAbsorptionCof[c1][0] == fSandiaTable[k2][0])
        {
	  flag = 1;
 	  break;
        }
      }
      if(flag == 0)
      {
        fPhotoAbsorptionCof[c][0] = fSandiaTable[k2][0];
	c++;
      }
    }       
  }   // end for(i)
  
  SandiaSort(fPhotoAbsorptionCof,c);
  fMaxInterval = c;
 
  const G4double* fractionW = fMaterial->GetFractionVector();
   
  for(i = 0; i < fMaxInterval; i++)
  {
      for(j = 1; j < 5; j++) fPhotoAbsorptionCof[i][j] = 0.;
  }
  for(i = 0; i < noElm; i++)
  {
      n1 = 1;
      G4double I1 = fIonizationPotentials[Z[i]]*keV;

      for(j = 1; j < Z[i]; j++)
      {
         n1 += fNbOfIntervals[j];
      }
      G4int n2 = n1 + fNbOfIntervals[Z[i]] - 1;

      for(k = n1; k < n2; k++)
      {
         G4double B1 = fSandiaTable[k][0];
         G4double B2 = fSandiaTable[k+1][0];
         for(G4int c = 1; c < fMaxInterval-1; c++)
         {
            G4double E1 = fPhotoAbsorptionCof[c][0];
            G4double E2 = fPhotoAbsorptionCof[c+1][0];
            if(B1 > E1 || B2 < E2 || E1 < I1)
	    {
	       continue;
	    }
	    for(j = 1; j < 5; j++)
  	    {
               fPhotoAbsorptionCof[c][j] += fSandiaTable[k][j]*fractionW[i];
 	    }
	  }  
       }   
       for(j = 1; j < 5; j++)   // Last interval
       {
          fPhotoAbsorptionCof[fMaxInterval-1][j] += fSandiaTable[k][j]*fractionW[i];
       }
  }     // for(i)

  c = 0;     // Deleting of first intervals where all coefficients = 0

  do                        
  {
    c++;

    if( fPhotoAbsorptionCof[c][1] != 0.0 ||
	fPhotoAbsorptionCof[c][2] != 0.0 ||
	fPhotoAbsorptionCof[c][3] != 0.0 || 
	fPhotoAbsorptionCof[c][4] != 0.0     )  continue;
       
    for(jj = 2; jj < fMaxInterval; jj++)
    {
      for(kk = 0; kk < 5; kk++)
      {
	     fPhotoAbsorptionCof[jj-1][kk]= fPhotoAbsorptionCof[jj][kk];
      }
    }
    fMaxInterval--;
    c--;
  }
  while(c < fMaxInterval - 1);
  	
  // create the sandia matrix for this material

  fMaxInterval--;  // vmg 20.11.10
    
  fMatSandiaMatrix = new G4OrderedTable();
 
  for (i = 0; i < fMaxInterval; i++)
  {
     fMatSandiaMatrix->push_back(new G4DataVector(5,0.));
  }	         	
  for ( i = 0; i < fMaxInterval; i++ )
  {
    for( j = 0; j < 5; j++ )
    {
      (*(*fMatSandiaMatrix)[i])[j] = fPhotoAbsorptionCof[i+1][j];
    }     
  }
  fMatNbOfIntervals = fMaxInterval; 
	         	    
  if ( fVerbose > 0 )
  {
    G4cout<<"vmg, G4SandiaTable::ComputeMatTable(), mat = "<<fMaterial->GetName()<<G4endl;

    for ( i = 0; i < fMaxInterval; i++ )
    {
      // G4cout<<i<<"\t"<<(*(*fMatSandiaMatrix)[i])[0]<<" keV \t"<<(*(*fMatSandiaMatrix)[i])[1]
      //       <<"\t"<<(*(*fMatSandiaMatrix)[i])[2]<<"\t"<<(*(*fMatSandiaMatrix)[i])[3]
      //   <<"\t"<<(*(*fMatSandiaMatrix)[i])[4]<<G4endl;

      G4cout<<i<<"\t"<<GetSandiaCofForMaterial(i,0)/keV<<" keV \t"<<this->GetSandiaCofForMaterial(i,1)
       <<"\t"<<this->GetSandiaCofForMaterial(i,2)<<"\t"<<this->GetSandiaCofForMaterial(i,3)
       <<"\t"<<this->GetSandiaCofForMaterial(i,4)<<G4endl;
    }
  }	         	    
  delete [] Z;
  return;
}  

//     G4SandiaTable class -- end of implementation file
//
////////////////////////////////////////////////////////////////////////////


