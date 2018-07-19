#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "PoPs.h"
#include "PoPs_mass.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static struct ZLabels {
    int Z;
    char const *Symbol;
} Zs[] = {  {   0,  "n" }, {   1,  "H" }, {   2, "He" }, {   3, "Li" }, {   4, "Be" }, {   5,  "B" }, {   6,  "C" }, {   7,  "N" }, {   8,  "O" },
            {   9,  "F" }, {  10, "Ne" }, {  11, "Na" }, {  12, "Mg" }, {  13, "Al" }, {  14, "Si" }, {  15,  "P" }, {  16,  "S" }, {  17, "Cl" },
            {  18, "Ar" }, {  19,  "K" }, {  20, "Ca" }, {  21, "Sc" }, {  22, "Ti" }, {  23,  "V" }, {  24, "Cr" }, {  25, "Mn" }, {  26, "Fe" },
            {  27, "Co" }, {  28, "Ni" }, {  29, "Cu" }, {  30, "Zn" }, {  31, "Ga" }, {  32, "Ge" }, {  33, "As" }, {  34, "Se" }, {  35, "Br" },
            {  36, "Kr" }, {  37, "Rb" }, {  38, "Sr" }, {  39,  "Y" }, {  40, "Zr" }, {  41, "Nb" }, {  42, "Mo" }, {  43, "Tc" }, {  44, "Ru" },
            {  45, "Rh" }, {  46, "Pd" }, {  47, "Ag" }, {  48, "Cd" }, {  49, "In" }, {  50, "Sn" }, {  51, "Sb" }, {  52, "Te" }, {  53,  "I" },
            {  54, "Xe" }, {  55, "Cs" }, {  56, "Ba" }, {  57, "La" }, {  58, "Ce" }, {  59, "Pr" }, {  60, "Nd" }, {  61, "Pm" }, {  62, "Sm" },
            {  63, "Eu" }, {  64, "Gd" }, {  65, "Tb" }, {  66, "Dy" }, {  67, "Ho" }, {  68, "Er" }, {  69, "Tm" }, {  70, "Yb" }, {  71, "Lu" },
            {  72, "Hf" }, {  73, "Ta" }, {  74,  "W" }, {  75, "Re" }, {  76, "Os" }, {  77, "Ir" }, {  78, "Pt" }, {  79, "Au" }, {  80, "Hg" },
            {  81, "Tl" }, {  82, "Pb" }, {  83, "Bi" }, {  84, "Po" }, {  85, "At" }, {  86, "Rn" }, {  87, "Fr" }, {  88, "Ra" }, {  89, "Ac" },
            {  90, "Th" }, {  91, "Pa" }, {  92,  "U" }, {  93, "Np" }, {  94, "Pu" }, {  95, "Am" }, {  96, "Cm" }, {  97, "Bk" }, {  98, "Cf" },
            {  99, "Es" }, { 100, "Fm" }, { 101, "Md" }, { 102, "No" }, { 103, "Lr" }, { 104, "Rf" }, { 105, "Db" }, { 106, "Sg" }, { 107, "Bh" },
            { 108, "Hs" }, { 109, "Mt" } };
static const int nZs = sizeof( Zs ) / sizeof( Zs[0] );

static char const *lPoPs_ZSymbol( int Z );
/*
========================================================================
*/
int lPoPs_addParticleIfNeeded( statusMessageReporting *smr, char const *name, char const *special ) {

    int index = PoPs_particleIndex( name ), ZA, Z = 0, A = 0,/* level = 0,*/ ispecial;
    char *endptr, name_[256], AStr[32];
    char const *ZStr, *alias = NULL;
    PoP *pop, *pop_;
 /* enum PoPs_genre genre = PoPs_genre_unknown; */
    char const *yiNames[]   = { "p",  "h2", "h3", "he3", "he4", "photon" };
    char const *yiAliases[] = { "h1",  "d",  "t", "he3",   "a", "g" };
 /* enum PoPs_genre yiGenres[] = { PoPs_genre_baryon, PoPs_genre_nucleus, PoPs_genre_nucleus, PoPs_genre_nucleus,
        PoPs_genre_nucleus, PoPs_genre_photon }; */

    if( special == NULL ) special = "";
    if( index < 0 ) {
        if( isdigit( name[0] ) ) {
            ZA = (int) strtol( name, &endptr, 10 );
            if( *endptr != 0 ) {
                smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "string '%s' not a value ZA", name );
                return( -1 );
            }
            Z = ZA / 1000;
            A = ZA % 1000;
            /*level = 0;*/
            ispecial = 0;
            if( strcmp( special, "LLNL" ) == 0 ) {
                if( ( ZA > 1 ) && ( ZA < 8 ) ) {
                    strcpy( name_, yiNames[ZA-2] );
                    alias = yiAliases[ZA-2];
                 /* genre = yiGenres[ZA-2];*/
                    ispecial = 1; }
                else if( ( ZA == 1801 ) || ( ZA == 1901 ) ) {
                    strcpy( name_, yiNames[0] );
                    alias = yiAliases[0];
                 /* genre = yiGenres[0]; */
                    ispecial = 1; }
                else if( ZA == 1902 ) {
                    strcpy( name_, yiNames[1] );
                    alias = yiAliases[1];
                 /* genre = yiGenres[1]; */
                    ispecial = 1; }
                else if( ZA == 4809 ) {
                    strcpy( name_, "Be9" );
                 /* genre = PoPs_genre_atom; */
                    ispecial = 1; }
                else if( ZA == 4909 ) {
                    strcpy( name_, "Be9" );
                 /* genre = PoPs_genre_atom; */
                    ispecial = 1; }
                else if( ZA == 6912 ) {
                    strcpy( name_, "C12" );
                 /* genre = PoPs_genre_atom; */
                    ispecial = 1; }
                else if( ZA == 8916 ) {
                    strcpy( name_, "O16" );
                 /* genre = PoPs_genre_atom; */
                    ispecial = 1; }
                else if( ZA == 95242 ) {
                    strcpy( name_, "Am242_e2" );
                    /*level = 2;*/
                 /* genre = PoPs_genre_atom; */
                    ispecial = 1; }
                else if( Z == 99 ) {
                    if( ( 120 <= A ) && ( A < 126 ) ) {
                        sprintf( name_, "FissionProductENDL99%d", A );
                     /* genre = PoPs_genre_atom; */
                        ispecial = 1;
                    }
                }
            }
            if( ispecial == 0 ) {
                if( ZA == 1 ) {
                    AStr[0] = 0; }
                else if( A == 0 ) {
                    strcpy( AStr, "_natural" ); }
                else {
                    sprintf( AStr, "%d", A );
                }
                if( ( ZStr = lPoPs_ZSymbol( Z ) ) == NULL ) {
                    smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "string '%s' not a value ZA; Z = %d is not supported", name, Z );
                    return( -1 );
                }
                sprintf( name_, "%s%s", ZStr, AStr );
              /* genre = PoPs_genre_atom; */
              /* if( ZA == 1 ) genre = PoPs_genre_baryon; */
            } }
        else {
            strcpy( name_, name );
            ZA = -1;
            if( strcmp( name, "neutron" ) == 0 ) {
                strcpy( name_, "n" );
                alias = name;
            /*  genre = PoPs_genre_baryon; */ }
            else if( strcmp( name, "electron" ) == 0 ) {
                strcpy( name_, "e-" );
                alias = name;
            /*  genre = PoPs_genre_lepton; */ } 
            else if( strcmp( name, "positron" ) == 0 ) {
                strcpy( name_, "e+" );
                alias = name;
            /*  genre = PoPs_genre_lepton; */ } 
            else if( ( strcmp( name, "h1" ) == 0 ) || ( strcmp( name, "proton" ) == 0 ) ) {
                ZA = 2; }
            else if( ( strcmp( name, "d" ) == 0 ) || ( strcmp( name, "deuteron" ) == 0 ) ) {
                ZA = 3; }
            else if( ( strcmp( name, "t" ) == 0 ) || ( strcmp( name, "triton" ) == 0 ) ) {
                ZA = 4; }
            else if( strcmp( name, "helium3" ) == 0 ) {
                ZA = 5; }
            else if( ( strcmp( name, "a" ) == 0 ) || ( strcmp( name, "alpha" ) == 0 ) || ( strcmp( name, "helium4" ) == 0 ) ) {
                ZA = 6; }
            else if( ( strcmp( name, "g" ) == 0 ) || ( strcmp( name, "gamma" ) == 0 ) ) {
                ZA = 7; }
            else if( strcmp( name, "FP" ) == 0 ) {
                        strcpy( name_, "FissionProductENDL99120" );
                    /*  genre = PoPs_genre_atom; */
            }
            if( ZA != -1 ) {
                strcpy( name_, yiNames[ZA-2] );
                alias = name;
             /* genre = yiGenres[ZA-2]; */
            }
        }

        if( ( pop = PoPs_particleCreateLoadInfo( smr, name_ ) ) == NULL ) {
            smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "particle '%s' converted to name '%s' not in database", name, name_ );
            return( -1 );
        }
        if( ( pop_ = PoPs_addParticleIfNeeded( smr, pop ) ) != pop ) PoP_free( pop );       /* Still need to add alias as index was < 0. */
        index = pop_->index;

        if( PoPs_particleIndex( name ) < 0 ) {
            if( ( pop = PoP_makeAlias( smr, name_, name ) ) == NULL ) return( -1 );
            if( ( pop_ = PoPs_addParticleIfNeeded( smr, pop ) ) != pop ) return( -1 );          /* pop_ should be pop as index was < 0. */
        }

        if( alias != NULL ) {
            if( PoPs_particleIndex( alias ) < 0 ) {
                if( ( pop = PoP_makeAlias( smr, name_, alias ) ) == NULL ) return( -1 );
                if( ( pop_ = PoPs_addParticleIfNeeded( smr, pop ) ) != pop ) return( -1 );  /* Required for some yis. */
            }
        }
    }
    return( index );
}
/*
========================================================================
*/
static char const *lPoPs_ZSymbol( int Z ) {

    //Coverity #63066
    if( ( Z < 0 ) || ( Z >= nZs ) ) return( NULL );
    return( Zs[Z].Symbol );
}

#if defined __cplusplus
}
#endif
