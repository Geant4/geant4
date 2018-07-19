// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_mpi_dummy_mpi_h
#define tools_mpi_dummy_mpi_h

extern "C" {

///////////////////////////////////////////////////////////////
/// to pass hd2mpi ////////////////////////////////////////////
///////////////////////////////////////////////////////////////

typedef void* MPI_Comm;
typedef void* MPI_Datatype;

#define MPI_UNSIGNED           0
#define MPI_FLOAT              0
#define MPI_DOUBLE             0
#define MPI_UNSIGNED_CHAR      0
#define MPI_CHAR               0
#define MPI_LONG               0
#define MPI_SHORT              0
#define MPI_INT                0
#define MPI_UNSIGNED_LONG      0
#define MPI_LONG_LONG          0
#define MPI_UNSIGNED_LONG_LONG 0

#define MPI_ANY_SOURCE 0

#define MPI_SUCCESS 1

#ifdef TOOLS_USE_MPI_PACK_NOT_CONST
inline int MPI_Pack(void*,int,MPI_Datatype,void*,int,int*,MPI_Comm){return 0;}
inline int MPI_Unpack(void*,int,int*,void*,int,MPI_Datatype,MPI_Comm){return 0;}
#else
inline int MPI_Pack(const void*,int,MPI_Datatype,void*,int,int*,MPI_Comm){return 0;}
inline int MPI_Unpack(const void*,int,int*,void*,int,MPI_Datatype,MPI_Comm){return 0;}
#endif

///////////////////////////////////////////////////////////////
/// to pass h2mpi, hs2mpi /////////////////////////////////////
///////////////////////////////////////////////////////////////

struct _MPI_Status {
  int MPI_SOURCE;
  int MPI_TAG;
};
typedef _MPI_Status MPI_Status;

inline int MPI_Probe(int,int,MPI_Comm,MPI_Status*){return 0;}
inline int MPI_Get_count(const MPI_Status*,MPI_Datatype,int*){return 0;}
inline int MPI_Send(const void*,int,MPI_Datatype,int,int,MPI_Comm){return 0;}
inline int MPI_Recv(void*,int,MPI_Datatype,int,int,MPI_Comm,MPI_Status*){return 0;}

///////////////////////////////////////////////////////////////
/// to pass examples/cpp/mpi.cpp //////////////////////////////
///////////////////////////////////////////////////////////////

#define MPI_COMM_WORLD 0
#define MPI_MAX_PROCESSOR_NAME 100

inline int MPI_Init(int*,char***){return 0;}
inline int MPI_Finalize(void){return 0;}
inline int MPI_Comm_size(MPI_Comm,int*){return 0;}
inline int MPI_Comm_rank(MPI_Comm,int*){return 0;}
inline int MPI_Get_processor_name(char*,int*){return 0;}

}


#endif
