// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SandiaTable.hh,v 1.7 2001-02-16 17:09:54 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// class description
//
// This class is an interface to G4StaticSandiaData.
// it provides - Sandia coeff for an element, given its Z
//             - sandia coeff for a material, given a pointer to it
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....
//
// 30.01.01 major bug in the computation of AoverAvo and in the units (/g!)
//          in GetSandiaCofPerAtom(). mma
// 18.11.98 simplified public interface; new methods for materials.  mma
// 10.06.97 created. V. Grichine
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

#ifndef G4SANDIATABLE_HH
#define G4SANDIATABLE_HH

#include "G4OrderedTable.hh"      
#include "G4ios.hh"
#include "globals.hh"
#include <assert.h>

class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

class G4SandiaTable
{
public:  // with description

    G4SandiaTable(G4int);	         
    G4SandiaTable(G4Material*);	         
   ~G4SandiaTable();

    //per atom of an Element:   
    static G4int     GetNbOfIntervals   (G4int Z);
    static G4double  GetSandiaCofPerAtom(G4int Z, G4int, G4int);
    static G4double* GetSandiaCofPerAtom(G4int Z, G4double energy);
    static G4double  GetIonizationPot   (G4int Z);
    static G4double  GetZtoA            (G4int Z);
    
    //per volume of a material:
           G4int     GetMatNbOfIntervals()  {return fMatNbOfIntervals;};
           G4double  GetSandiaCofForMaterial(G4int,G4int);
           G4double* GetSandiaCofForMaterial(G4double energy);
	   

private:
       
    void ComputeMatSandiaMatrix();
    
private:

    static const G4double        fSandiaTable[981][5];
    static const G4int           fNbOfIntervals[101];
    static       G4int           fCumulInterval[101];
    static const G4double        fZtoAratio[101];
    static const G4double        fIonizationPotentials[101];
               
    static       G4double        fSandiaCofPerAtom[4];
    
                 G4Material*     fMaterial;
                 G4int           fMatNbOfIntervals;
                 G4OrderedTable* fMatSandiaMatrix;
		   

/////////////////////////////////////////////////////////////////////
//
// Methods for PAI model

public:  // without description

         inline void SandiaSwap( G4double** da,
                                 G4int i,
                                 G4int j );

         void SandiaSort( G4double** da,
                          G4int sz );

	 G4int SandiaIntervals( G4int Z[],
                                G4int el );

         G4int SandiaMixing(       G4int Z[],
                             const G4double fractionW[],
                                   G4int el,
                                   G4int mi );

         inline G4double GetPhotoAbsorpCof(G4int i , G4int j)const ;

         inline G4int GetMaxInterval() const ;



//////////////////////////////////////////////////////////////////////////
//
// data members for PAI model

private:
         static const G4int    fNumberOfElements  ;
         static const G4int    fIntervalLimit ;
         static const G4int    fNumberOfIntervals  ;
		
         // G4double fPhotoAbsorptionCof[101][5] ; // SandiaTable  for mixture

	 G4double** fPhotoAbsorptionCof ;	// SandiaTable  for mixture

         // G4OrderedTable*
	 G4int fMaxInterval ;

//
//
//////////////////////////////////////////////////////////////////////////
  
};
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

// Inline function implementations


inline
G4int G4SandiaTable::GetNbOfIntervals(G4int Z)
{
   return fNbOfIntervals[Z];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4double  G4SandiaTable::GetSandiaCofPerAtom(G4int Z, G4int interval, G4int j)
{
   assert (Z>0 && Z<101 && interval>=0 && interval<fNbOfIntervals[Z]
                        && j>=0 && j<5);

   G4int row = fCumulInterval[Z-1] + interval;
   if (j==0) return fSandiaTable[row][0]*keV;
   G4double AoverAvo = Z*amu/fZtoAratio[Z];         
   return AoverAvo*(fSandiaTable[row][j]*cm2*pow(keV,G4double(j))/g);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4double* G4SandiaTable::GetSandiaCofPerAtom(G4int Z, G4double energy)
{
   assert (Z > 0);
   
   G4double Emin  = fSandiaTable[fCumulInterval[Z-1]][0]*keV;
   G4double Iopot = fIonizationPotentials[Z]*eV;
   if (Iopot > Emin) Emin = Iopot;
   
   G4int interval = fNbOfIntervals[Z] - 1;
   G4int row = fCumulInterval[Z-1] + interval;
   while ((interval>0) && (energy<fSandiaTable[row][0]*keV))
         row = fCumulInterval[Z-1] + --interval;

   if (energy >= Emin)
     {        
      G4double AoverAvo = Z*amu/fZtoAratio[Z];
         
      fSandiaCofPerAtom[0]=AoverAvo*(fSandiaTable[row][1]*cm2*keV/g);     
      fSandiaCofPerAtom[1]=AoverAvo*(fSandiaTable[row][2]*cm2*keV*keV/g);     
      fSandiaCofPerAtom[2]=AoverAvo*(fSandiaTable[row][3]*cm2*keV*keV*keV/g);     
      fSandiaCofPerAtom[3]=AoverAvo*(fSandiaTable[row][4]*cm2*keV*keV*keV*keV/g);
     }
   else
      fSandiaCofPerAtom[0] = fSandiaCofPerAtom[1] = fSandiaCofPerAtom[2] =
      fSandiaCofPerAtom[3] = 0.;
                        
   return fSandiaCofPerAtom;     

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4double G4SandiaTable::GetIonizationPot(G4int Z)
{
   return fIonizationPotentials[Z]*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

inline
G4double G4SandiaTable::GetZtoA(G4int Z)
{
   return fZtoAratio[Z];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... ....oooOO0OOooo....

///////////////////////////////////////////////////////////////////////
//
// Inline methods for PAI model

inline
void
G4SandiaTable::SandiaSwap( G4double** da ,
                           G4int i,
                           G4int j )
{
  G4double tmp = da[i][0] ;
  da[i][0] = da[j][0] ;
  da[j][0] = tmp ;
}

/////////////////////////////////////////////////////////////////////////
//
//

inline
G4double G4SandiaTable::GetPhotoAbsorpCof(G4int i, G4int j) const
{
   G4double unitCof ;
   if(j == 0)
   {
      unitCof = keV ;
   }
   else
   {
      unitCof = (cm2/g)*pow(keV,(G4double)j) ;
   }
   return  fPhotoAbsorptionCof[i][j]*unitCof ;
}

/////////////////////////////////////////////////////////////////////
//
//

inline
G4int G4SandiaTable::GetMaxInterval() const
{
   return fMaxInterval ;
}

//
//
////////////////////////////////////////////////////////////////////////////


#endif 

