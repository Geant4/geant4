/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef PoPs_h_included
#define PoPs_h_included

/* Disable Effective C++ warnings in PoP code. */
#if __INTEL_COMPILER > 1399
#pragma warning( disable:593 )
#endif

#include <statusMessageReporting.h>
/*
* MPI stuff.
*/
#ifdef PoPs_MPI
#include <mpi.h>
#endif

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

#define POPS_VERSION_MAJOR 1
#define POPS_VERSION_MINOR 0
#define POPS_VERSION_PATCHLEVEL 5

#define PoPs_packageSymbol "PoPs (properties of particles)"
#define PoPs_packageName PoPs_packageSymbol " (properties of particles)"
typedef struct PoP_s PoP;

enum PoPs_errorTokens { PoPs_errorToken_Okay, PoPs_errorToken_badName, PoPs_errorToken_badIndex, PoPs_errorToken_badUnitConversion };
enum PoPs_genre { PoPs_genre_invalid, PoPs_genre_unknown, PoPs_genre_alias, PoPs_genre_photon, PoPs_genre_lepton, 
        PoPs_genre_quark, PoPs_genre_meson, PoPs_genre_baryon, PoPs_genre_nucleus, PoPs_genre_atom };
/*
*   In the following struct, 'index' is the index of the particle (proper or aliased) in the list of particles.  If a particle 
*   is a proper particle its properIndex is -1. Otherwise, it is the index of the aliased particle's proper particle. If a proper 
*   particle does not have an aliased particle referring to it, aliasIndex is -1. If a proper particle has aliaes particles, 
*   its aliasIndex is the index of its first aliased particle. If a second alias is added to a proper particle, then its first 
*   aliased particle's aliasIndex is the index of that particle, and so on. The last aliased particle added has aliasIndex = -1.
*/
struct PoP_s {              /* Any changes here must be reflected in functions PoP_initialize and PoP_copyParticle and in file PoPs_Bcast.c logic. */
    int index, properIndex, aliasIndex;
    enum PoPs_genre genre;
    char const *name;
    int Z, A, l;
    double mass;            /* Mass to be added to base. */
    char const *massUnit;
};

extern int PoPs_smr_ID;

const char *PoPs_version( void );
int PoPs_versionMajor( void );
int PoPs_versionMinor( void );
int PoPs_versionPatchLevel( void );

int PoPs_register( void );
int PoPs_readDatabase( statusMessageReporting *smr, char const *fileName );
int PoPs_release( statusMessageReporting *smr );
PoP *PoPs_addParticleIfNeeded( statusMessageReporting *smr, PoP *pop );
PoP *PoPs_copyAddParticleIfNeeded( statusMessageReporting *smr, PoP *pop );
PoP *PoPs_addAliasIfNeeded( statusMessageReporting *smr, char const *name, char const *alias );
int PoPs_numberOfParticle( void );
int PoPs_particleIndex( char const *name );
int PoPs_particleIndex_smr( statusMessageReporting *smr, char const *name, char const *file, int line, char const *func );
char const *PoPs_getName_atIndex( statusMessageReporting *smr, int index );
double PoPs_getMassInUnitOf( statusMessageReporting *smr, char const *name, char const *unit );
double PoPs_getMassInUnitOf_atIndex( statusMessageReporting *smr, int index, char const *unit );
enum PoPs_genre PoPs_getGenre( statusMessageReporting *smr, char const *name );
enum PoPs_genre PoPs_getGenre_atIndex( statusMessageReporting *smr, int index );
int PoPs_getZ_A_l( statusMessageReporting *smr, char const *name, int *Z, int *A, int *l );
int PoPs_getZ_A_l_atIndex( statusMessageReporting *smr, int index, int *Z, int *A, int *l );
int PoPs_hasNucleus( statusMessageReporting *smr, char const *name, int protonIsNucleus );
int PoPs_hasNucleus_atIndex( statusMessageReporting *smr, int index, int protonIsNucleus );
char const *PoPs_getAtomsName( statusMessageReporting *smr, char const *name );
char const *PoPs_getAtomsName_atIndex( statusMessageReporting *smr, int index );
int PoPs_getAtomsIndex( statusMessageReporting *smr, char const *name );
int PoPs_getAtomsIndex_atIndex( statusMessageReporting *smr, int index );
PoP *PoPs_getParticle_atIndex( int index );

char const *PoPs_genreTokenToString( enum PoPs_genre genre );
void PoPs_print( int sorted );
void PoPs_write( FILE *f, int sorted );

PoP *PoP_new( statusMessageReporting *smr );
int PoP_initialize( statusMessageReporting *smr, PoP *pop );
int PoP_release( PoP *pop );
PoP *PoP_free( PoP *pop );
int PoP_copyParticle( statusMessageReporting *smr, PoP *desc, PoP *src );
PoP *PoP_makeParticle( statusMessageReporting *smr, enum PoPs_genre genre, char const *name, double mass, char const *massUnit );
int PoP_setZ_A_l( statusMessageReporting *smr, PoP *pop, int Z, int A, int l );
int PoP_getIndex( PoP *pop );
char const *PoP_getName( PoP *pop );

int PoPs_particleReadDatabase( statusMessageReporting *smr, char const *name );
PoP *PoPs_particleCreateLoadInfo( statusMessageReporting *smr, const char *name );
int PoPs_particleLoadInfo( statusMessageReporting *smr, const char *name, PoP *pop );

double PoP_getMassInUnitOf( statusMessageReporting *smr, PoP *pop, char const *unit );

PoP *PoP_makeAlias( statusMessageReporting *smr, char const *name, char const *alias );

int PoPs_unitConversionRatio( char const *_from, char const *_to, double *ratio );

int lPoPs_addParticleIfNeeded( statusMessageReporting *smr, char const *name, char const *special );

/*
* MPI stuff.
*/
#ifdef PoPs_MPI
int PoPs_Bcast( statusMessageReporting *smr, MPI_Comm comm, int bossRank );
#endif

/* Use the next function with caution as it is only for initial testing of the package and will soon be gone. */
int PoPs_setBDFLS_File( char const *name );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of PoPs_h_included. */
