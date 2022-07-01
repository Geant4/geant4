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
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 21-04-16, created by E.Bagli

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4CrystalUnitCell.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>

G4CrystalUnitCell::G4CrystalUnitCell(G4double sizeA,
                                     G4double sizeB,
                                     G4double sizeC,
                                     G4double alpha,
                                     G4double beta,
                                     G4double gamma,
                                     G4int spacegroup):
theSpaceGroup(spacegroup),
theSize(G4ThreeVector(sizeA,sizeB,sizeC)),
theAngle(G4ThreeVector(alpha,beta,gamma))
{
    
    nullVec = G4ThreeVector(0.,0.,0.);
    theUnitBasis[0] = CLHEP::HepXHat;
    theUnitBasis[1] = CLHEP::HepYHat;
    theUnitBasis[2] = CLHEP::HepZHat;

    theRecUnitBasis[0] = CLHEP::HepXHat;
    theRecUnitBasis[1] = CLHEP::HepYHat;
    theRecUnitBasis[2] = CLHEP::HepZHat;

    cosa=std::cos(alpha), cosb=std::cos(beta), cosg=std::cos(gamma);
    sina=std::sin(alpha), sinb=std::sin(beta), sing=std::sin(gamma);

    cosar = (cosb*cosg-cosa)/(sinb*sing);
    cosbr = (cosa*cosg-cosb)/(sina*sing);
    cosgr = (cosa*cosb-cosg)/(sina*sinb);
    
    theVolume = ComputeCellVolume();
    theRecVolume = 1. / theVolume;

    theRecSize[0] = sizeB * sizeC * sina / theVolume;
    theRecSize[1] = sizeC * sizeA * sinb / theVolume;
    theRecSize[2] = sizeA * sizeB * sing / theVolume;
    
    theRecAngle[0] = std::acos(cosar);
    theRecAngle[1] = std::acos(cosbr);
    theRecAngle[2] = std::acos(cosgr);
    
    G4double x3,y3,z3;
    
    switch (GetLatticeSystem(theSpaceGroup)) {
        case Amorphous:
            break;
        case Cubic: // Cubic, C44 set
            break;
        case Tetragonal:
            break;
        case Orthorhombic:
            break;
        case Rhombohedral:
            theUnitBasis[1].rotateZ(gamma-CLHEP::halfpi);	// X-Y opening angle
            // Z' axis computed by hand to get both opening angles right
            // X'.Z' = cos(alpha), Y'.Z' = cos(beta), solve for Z' components
            x3=cosa, y3=(cosb-cosa*cosg)/sing, z3=std::sqrt(1.-x3*x3-y3*y3);
            theUnitBasis[2] = G4ThreeVector(x3, y3, z3).unit();
            break;
        case Monoclinic:
            theUnitBasis[2].rotateX(beta-CLHEP::halfpi);	// Z-Y opening angle
            break;
        case Triclinic:
            theUnitBasis[1].rotateZ(gamma-CLHEP::halfpi);	// X-Y opening angle
            // Z' axis computed by hand to get both opening angles right
            // X'.Z' = cos(alpha), Y'.Z' = cos(beta), solve for Z' components
            x3=cosa, y3=(cosb-cosa*cosg)/sing, z3=std::sqrt(1.-x3*x3-y3*y3);
            theUnitBasis[2] = G4ThreeVector(x3, y3, z3).unit();
            break;
        case Hexagonal:  // Tetragonal, C16=0
            theUnitBasis[1].rotateZ(30.*CLHEP::deg);	// X-Y opening angle
            break;
        default:
            break;
    }
    
    for(auto i:{0,1,2}){
        theBasis[i] = theUnitBasis[i] * theSize[i];
        theRecBasis[i] = theRecUnitBasis[i] * theRecSize[i];
    }
    
    // Initialize sgInfo
    /* at first some initialization for SgInfo */
    /*
    const T_TabSgName *tsgn = NULL;
    
    SgInfo.MaxList = 192;
    SgInfo.ListSeitzMx = malloc( SgInfo.MaxList * sizeof(*SgInfo.ListSeitzMx) );

    // no list info needed here
    SgInfo.ListRotMxInfo = NULL;
    tsgn = FindTabSgNameEntry(SchoenfliesSymbols[theSpaceGroup], 'A');

    // initialize SgInfo struct
    InitSgInfo( &SgInfo );
    SgInfo.TabSgName = tsgn;
    if ( tsgn ){
        SgInfo.GenOption = 1;
    }
    
    ParseHallSymbol( SchoenfliesSymbols[theSpaceGroup], &SgInfo );
    CompleteSgInfo( &SgInfo );
    Set_si( &SgInfo );
    */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CrystalUnitCell::~G4CrystalUnitCell(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

theLatticeSystemType G4CrystalUnitCell::GetLatticeSystem(G4int aGroup){
    
    if(     aGroup >=   1 && aGroup <=  2 )  {return Triclinic;}
    else if(aGroup >=   3 && aGroup <= 15 )  {return Monoclinic;}
    else if(aGroup >=  16 && aGroup <= 74 )  {return Orthorhombic;}
    else if(aGroup >=  75 && aGroup <= 142)  {return Tetragonal;}
    else if(aGroup == 146 || aGroup == 148 ||
            aGroup == 155 || aGroup == 160 ||
            aGroup == 161 || aGroup == 166 ||
            aGroup == 167)                   {return Rhombohedral;}
    else if(aGroup >= 143 && aGroup <= 167)  {return Hexagonal;}
    else if(aGroup >= 168 && aGroup <= 194)  {return Hexagonal;}
    else if(aGroup >= 195 && aGroup <= 230)  {return Cubic;}
    
    return Amorphous;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
theBravaisLatticeType G4CrystalUnitCell::GetBravaisLattice(G4int aGroup){
    ;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ThreeVector& G4CrystalUnitCell::GetUnitBasis(G4int idx) const {
    return (idx>=0 && idx<3 ? theUnitBasis[idx] : nullVec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ThreeVector& G4CrystalUnitCell::GetBasis(G4int idx) const {
    return (idx>=0 && idx<3 ? theBasis[idx] : nullVec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ThreeVector& G4CrystalUnitCell::GetRecUnitBasis(G4int idx) const {
    return (idx>=0 && idx<3 ? theRecUnitBasis[idx] : nullVec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ThreeVector& G4CrystalUnitCell::GetRecBasis(G4int idx) const {
    return (idx>=0 && idx<3 ? theRecBasis[idx] : nullVec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector G4CrystalUnitCell::GetUnitBasisTrigonal(){
    // Z' axis computed by hand to get both opening angles right
    // X'.Z' = cos(alpha), Y'.Z' = cos(beta), solve for Z' components
    G4double x3=cosa, y3=(cosb-cosa*cosg)/sing, z3=std::sqrt(1.-x3*x3-y3*y3);
    return G4ThreeVector(x3, y3, z3).unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillAtomicUnitPos(G4ThreeVector& pos, std::vector<G4ThreeVector>& vecout){
    // Just for testing the infrastructure
    G4ThreeVector aaa = pos;
    vecout.push_back(aaa);
    vecout.emplace_back(2., 5., 3.);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillAtomicPos(G4ThreeVector& posin, std::vector<G4ThreeVector>& vecout){
    FillAtomicUnitPos(posin,vecout);
    for(auto &vec:vecout){
        vec.setX(vec.x()*theSize[0]);
        vec.setY(vec.y()*theSize[1]);
        vec.setZ(vec.z()*theSize[2]);
    }
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillElReduced(G4double Cij[6][6]) {
    switch (GetLatticeSystem()) {
        case Amorphous:
            return FillAmorphous(Cij);
            break;
        case Cubic: // Cubic, C44 set
            return FillCubic(Cij);
            break;
        case Tetragonal:
            return FillTetragonal(Cij);
            break;
        case Orthorhombic:
            return FillOrthorhombic(Cij);
            break;
        case Rhombohedral:
            return FillRhombohedral(Cij);
            break;
        case Monoclinic:
            return FillMonoclinic(Cij);
            break;
        case Triclinic:
            return FillTriclinic(Cij);
            break;
        case Hexagonal:  // Tetragonal, C16=0
            return FillHexagonal(Cij);
            break;
        default:
            break;
    }
    
    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillAmorphous(G4double Cij[6][6]) const {
    Cij[3][3] = 0.5*(Cij[0][0]-Cij[0][1]);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillCubic(G4double Cij[6][6]) const {
    G4double C11=Cij[0][0], C12=Cij[0][1], C44=Cij[3][3];
    
    for (size_t i=0; i<6; i++) {
        for (size_t j=i; j<6; j++) {
          if(i < 3 && j < 3)
          {
            Cij[i][j] = (i == j) ? C11 : C12;
          }
          else if(i == j && i >= 3)
          {
            Cij[i][i] = C44;
          }
          else
          {
            Cij[i][j] = 0.;
          }
        }
    }
    
    ReflectElReduced(Cij);
    
    return (C11!=0. && C12!=0. && C44!=0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillTetragonal(G4double Cij[6][6]) const {
    G4double C11=Cij[0][0], C12=Cij[0][1], C13=Cij[0][2], C16=Cij[0][5];
    G4double C33=Cij[2][2], C44=Cij[3][3], C66=Cij[5][5];
    
    Cij[1][1] = C11;	// Copy small number of individual elements
    Cij[1][2] = C13;
    Cij[1][5] = -C16;
    Cij[4][4] = C44;
    
    ReflectElReduced(Cij);
    
    // NOTE:  Do not test for C16 != 0., to allow calling from Hexagonal
    return (C11!=0. && C12!=0. && C13!=0. && C33!=0. && C44!=0. && C66!=0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillOrthorhombic(G4double Cij[6][6]) const {
    // No degenerate elements; just check for all non-zero
    ReflectElReduced(Cij);
    
    G4bool good = true;
    for (size_t i=0; i<6; i++) {
      for(size_t j = i + 1; j < 3; j++)
      {
        good &= (Cij[i][j] != 0);
      }
    }
    
    return good;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillRhombohedral(G4double Cij[6][6]) const {
    G4double C11=Cij[0][0], C12=Cij[0][1], C13=Cij[0][2], C14=Cij[0][3];
    G4double C15=Cij[0][4], C33=Cij[2][2], C44=Cij[3][3], C66=0.5*(C11-C12);
    
    Cij[1][1] = C11;	// Copy small number of individual elements
    Cij[1][2] = C13;
    Cij[1][3] = -C14;
    Cij[1][4] = -C15;
    Cij[3][5] = -C15;
    Cij[4][4] = C44;
    Cij[4][5] = C14;
    
    // NOTE:  C15 may be zero (c.f. rhombohedral(I) vs. (II))
    return (C11!=0 && C12!=0 && C13!=0 && C14!=0. &&
            C33!=0. && C44!=0. && C66!=0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillMonoclinic(G4double Cij[6][6]) const {
    // The monoclinic matrix has 13 independent elements with no degeneracies
    // Sanity condition is same as orthorhombic, plus C45, C(1,2,3)6
    
    return (FillOrthorhombic(Cij) && Cij[0][5]!=0. && Cij[1][5]!=0. &&
            Cij[2][5] != 0. && Cij[3][4]!=0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillTriclinic(G4double Cij[6][6]) const {
    // The triclinic matrix has the entire upper half filled (21 elements)
    
    ReflectElReduced(Cij);
    
    G4bool good = true;
    for (size_t i=0; i<6; i++) {
      for(size_t j = i; j < 6; j++)
      {
        good &= (Cij[i][j] != 0);
      }
    }
    
    return good;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::FillHexagonal(G4double Cij[6][6]) const {
    Cij[0][5] = 0.;
    Cij[4][5] = 0.5*(Cij[0][0] - Cij[0][1]);
    return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CrystalUnitCell::ReflectElReduced(G4double Cij[6][6]) const {
    for (size_t i=1; i<6; i++) {
        for (size_t j=i+1; j<6; j++) {
            Cij[j][i] = Cij[i][j];
        }
    }
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4CrystalUnitCell::ComputeCellVolume(){
    G4double a = theSize[0], b = theSize[1], c = theSize[2];
    
    switch(GetLatticeSystem())
    {
        case Amorphous:
            return 0.;
            break;
        case Cubic:
            return a * a * a;
            break;
        case Tetragonal:
            return a * a * c;
            break;
        case Orthorhombic:
            return a * b * c;
            break;
        case Rhombohedral:
            return a*a*a*std::sqrt(1.-3.*cosa*cosa+2.*cosa*cosa*cosa);
            break;
        case Monoclinic:
            return a*b*c*sinb;
            break;
        case Triclinic:
            return a*b*c*std::sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg*2.*cosa*cosb*cosg);
            break;
        case Hexagonal:
            return  std::sqrt(3.0)/2.*a*a*c;
            break;
        default:
            break;
    }
    
    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrystalUnitCell::GetIntSp2(G4int h,
                                      G4int k,
                                      G4int l){
    
    /* Reference:
     Table 2.4, pag. 65
     
     @Inbook{Ladd2003,
     author="Ladd, Mark and Palmer, Rex",
     title="Lattices and Space-Group Theory",
     bookTitle="Structure Determination by X-ray Crystallography",
     year="2003",
     publisher="Springer US",
     address="Boston, MA",
     pages="51--116",
     isbn="978-1-4615-0101-5",
     doi="10.1007/978-1-4615-0101-5_2",
     url="http://dx.doi.org/10.1007/978-1-4615-0101-5_2"
     }
    */
    
    G4double a = theSize[0], b = theSize[1], c = theSize[2];
    G4double a2 = a*a, b2 = b*b, c2 = c*c;
    G4double h2 = h*h, k2 = k*k, l2 = l*l;
    
    G4double cos2a,sin2a,sin2b;
    G4double R,T;
    
    switch(GetLatticeSystem())
    {
        case Amorphous:
            return 0.;
            break;
        case Cubic:
            return a2 / ( h2+k2+l2 );
            break;
        case Tetragonal:
            return 1.0 / ( (h2 + k2)/a2 + l2/c2 );
            break;
        case Orthorhombic:
            return 1.0 / ( h2/a2 + k2/b2 + l2/c2 );
            break;
        case Rhombohedral:
            cos2a=cosa*cosa; sin2a=sina*sina;
            T = h2+k2+l2+2.*(h*k+k*l+h*l) * ((cos2a-cosa)/sin2a);
            R = sin2a / (1. - 3*cos2a + 2.*cos2a*cosa);
            return a*a / (T*R);
            break;
        case Monoclinic:
            sin2b=sinb*sinb;
            return 1./(1./sin2b * (h2/a2+l2/c2-2*h*l*cosb/(a*c)) + k2/b2);
            break;
        case Triclinic:
            return 1./GetRecIntSp2(h,k,l);
            break;
        case Hexagonal:
            return 1. / ( (4.*(h2+k2+h*k) / (3.*a2)) + l2/c2 );
            break;
        default:
            break;
    }

    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrystalUnitCell::GetRecIntSp2(G4int h,
                                         G4int k,
                                         G4int l){
    /* Reference:
     Table 2.4, pag. 65
     
     @Inbook{Ladd2003,
     author="Ladd, Mark and Palmer, Rex",
     title="Lattices and Space-Group Theory",
     bookTitle="Structure Determination by X-ray Crystallography",
     year="2003",
     publisher="Springer US",
     address="Boston, MA",
     pages="51--116",
     isbn="978-1-4615-0101-5",
     doi="10.1007/978-1-4615-0101-5_2",
     url="http://dx.doi.org/10.1007/978-1-4615-0101-5_2"
     }
     */

    G4double a = theRecSize[0], b = theRecSize[1], c = theRecSize[2];
    G4double a2 = a*a, b2 = b*b, c2 = c*c;
    G4double h2 = h*h, k2 = k*k, l2 = l*l;

    switch(GetLatticeSystem())
    {
        case Amorphous:
            return 0.;
            break;
        case Cubic:
            return a2 * (h2+k2+l2);
            break;
        case Tetragonal:
            return (h2+k2)*a2 + l2*c2 ;
            break;
        case Orthorhombic:
            return h2*a2 + k2+b2 + h2*c2;
            break;
        case Rhombohedral:
            return (h2+k2+l2+2.*(h*k+k*l+h*l) * cosar)*a2;
            break;
        case Monoclinic:
            return h2*a2+k2*b2+l2*c2+2.*h*l*a*c*cosbr;
            break;
        case Triclinic:
            return h2*a2+k2*b2+l2*c2+2.*k*l*b*c*cosar+2.*l*h*c*a*cosbr+2.*h*k*a*b*cosgr;
            break;
        case Hexagonal:
            return (h2+k2+h*k)*a2 + l2*c2;
            break;
        default:
            break;
    }
    
    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrystalUnitCell::GetIntCosAng(G4int h1,
                                         G4int k1,
                                         G4int l1,
                                         G4int h2,
                                         G4int k2,
                                         G4int l2){

    /* Reference:
     Table 2.4, pag. 65
     
     @Inbook{Kelly2012,
     author="Anthony A. Kelly and Kevin M. Knowles",
     title="Appendix 3 Interplanar Spacings and Interplanar Angles",
     bookTitle="Crystallography and Crystal Defects, 2nd Edition",
     year="2012",
     publisher="John Wiley & Sons, Ltd.",
     isbn="978-0-470-75014-8",
     doi="10.1002/9781119961468",
     url="http://onlinelibrary.wiley.com/book/10.1002/9781119961468"
     }
     */
    
    G4double a = theRecSize[0], b = theRecSize[1], c = theRecSize[2];
    G4double a2 = a*a, b2 = b*b, c2 = c*c;
    G4double dsp1dsp2;
    switch(GetLatticeSystem())
    {
        case Amorphous:
            return 0.;
            break;
        case Cubic:
            return (h1*h2 + k1*k2 + l1+l2) / (std::sqrt(h1*h1 + k1*k1 + l1*l1) * std::sqrt(h2*h2 + k2*k2 + l2*l2));
            break;
        case Tetragonal:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return 0. ;
            break;
        case Orthorhombic:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return dsp1dsp2 * (h1*h2*a2 + k1*k2*a2 + l1*l2*c2);
            break;
        case Rhombohedral:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return dsp1dsp2 * (h1*h2*a2 + k1*k2*b2 + l1*l2*c2+
                               (k1*l2+k2*l1)*b*c*cosar+
                               (h1*l2+h2*l1)*a*c*cosbr+
                               (h1*k2+h2*k1)*a*b*cosgr);
            break;
        case Monoclinic:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return dsp1dsp2 * (h1*h2*a2 + k1*k2*b2 + l1*l2*c2+
                               (k1*l2+k2*l1)*b*c*cosar+
                               (h1*l2+h2*l1)*a*c*cosbr+
                               (h1*k2+h2*k1)*a*b*cosgr);
            break;
        case Triclinic:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return dsp1dsp2 * (h1*h2*a2 + k1*k2*b2 + l1*l2*c2+
                               (k1*l2+k2*l1)*b*c*cosar+
                               (h1*l2+h2*l1)*a*c*cosbr+
                               (h1*k2+h2*k1)*a*b*cosgr);
            break;
        case Hexagonal:
            dsp1dsp2 = std::sqrt(GetIntSp2(h1,k1,l1)*GetIntSp2(h2,k2,l2));
            return dsp1dsp2 *( (h1*h2 + k1*k2 + 0.5*(h1*k2+k1*h2))*a2 + l1*l2*c2);
            break;
        default:
            break;
    }

    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

