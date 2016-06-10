/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef PoPs_Bcast_private_h_included
#define PoPs_Bcast_private_h_included

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif


int PoPs_Bcast2( statusMessageReporting *smr, MPI_Comm comm, int bossRank, unitsDB *unitsRoot, PoPs *popsRoot );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of PoPs_Bcast_private_h_included. */
