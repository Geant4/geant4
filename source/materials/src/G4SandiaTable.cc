// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SandiaTable.cc,v 1.1 1999-01-07 16:09:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// 18.11.98 simplified public interface; new methods for materials.  mma
// 10.06.97 created. V. Grichine
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....


#include "G4SandiaTable.hh"
#include "G4StaticSandiaData.hh"
#include "G4Material.hh"

G4int    G4SandiaTable::fCumulInterval[101];
G4double G4SandiaTable::fSandiaCofPerAtom[4];
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4SandiaTable::G4SandiaTable(G4Material* material)
:fMaterial(material)
{
  //build the CumulInterval array
  fCumulInterval[0] = 1;
  for (G4int Z=1; Z<101; Z++)
      fCumulInterval[Z] = fCumulInterval[Z-1] + fNbOfIntervals[Z];
  
  //compute macroscopic Sandia coefs for a material   
  ComputeMatSandiaMatrix();    
}
							
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4SandiaTable::~G4SandiaTable()
{ delete fMatSandiaMatrix;}
						 	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

void G4SandiaTable::ComputeMatSandiaMatrix()
{  
  //get list of elements
  const G4int NbElm = fMaterial->GetNumberOfElements();
  const G4ElementVector* ElementVector = fMaterial->GetElementVector();
  
  G4int* Z = new G4int[NbElm];               //Atomic number
     
  //determine the total number of energy-intervals for this material
  fMatNbOfIntervals = 0;
  G4int elm;    
  for (elm=0; elm<NbElm; elm++)
     { Z[elm] = (int)(*ElementVector)(elm)->GetZ();
       fMatNbOfIntervals += fNbOfIntervals[Z[elm]];
     }  
     
  //create the sandia matrix for this material  
  fMatSandiaMatrix = new G4OrderedTable(fMatNbOfIntervals);
  G4int interval;
  for (interval=0; interval<fMatNbOfIntervals; interval++)
     (*fMatSandiaMatrix)(interval) = new G4ValVector(5);
  
  //copy the Energy bins (take care of the Ionization Potential)
  G4double Ebin;
  interval=0;
  for (elm=0; elm<NbElm; elm++)
     for (G4int row=fCumulInterval[Z[elm]-1];row<fCumulInterval[Z[elm]];row++)
        { Ebin = fSandiaTable[row][0]*keV;
          if ((row==fCumulInterval[Z[elm]-1])&&(GetIonizationPot(Z[elm])<Ebin))
              Ebin = GetIonizationPot(Z[elm]);
          (*(*fMatSandiaMatrix)(interval++))(0) = Ebin;
        }  
         
  //sort the energies in increasing values
  G4double tmp;
  for (G4int i1=0; i1<fMatNbOfIntervals; i1++)
     for (G4int i2=i1+1; i2<fMatNbOfIntervals; i2++)
        {if ((*(*fMatSandiaMatrix)(i1))(0) > (*(*fMatSandiaMatrix)(i2))(0))        
           {
            tmp = (*(*fMatSandiaMatrix)(i1))(0);
            (*(*fMatSandiaMatrix)(i1))(0) = (*(*fMatSandiaMatrix)(i2))(0);
            (*(*fMatSandiaMatrix)(i2))(0) = tmp;
           }
        }
        
  //ready to compute the Sandia coefs for the material
  const G4double* NbOfAtomsPerVolume = fMaterial->GetVecNbOfAtomsPerVolume();
  for (interval=0; interval<fMatNbOfIntervals; interval++)
     {
      Ebin = (*(*fMatSandiaMatrix)(interval))(0);        
      for (elm=0; elm<NbElm; elm++)
         {
           GetSandiaCofPerAtom(Z[elm], Ebin);
           for (G4int j=1; j<5; j++)
              (*(*fMatSandiaMatrix)(interval))(j) += NbOfAtomsPerVolume[elm]*
                                                     fSandiaCofPerAtom[j-1];
         }
     }
  delete [] Z;              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double  G4SandiaTable::GetSandiaCofForMaterial(G4int interval, G4int j)                                                 
{
   assert (interval>=0 && interval<fMatNbOfIntervals && j>=0 && j<5);                      
   return ((*(*fMatSandiaMatrix)(interval))(j)); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

G4double* G4SandiaTable::GetSandiaCofForMaterial(G4double energy)
{
   G4int interval = fMatNbOfIntervals - 1;
   while ((interval>0)&&(energy<(*(*fMatSandiaMatrix)(interval))(0))) interval--; 
   return &((*(*fMatSandiaMatrix)(interval))(1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
