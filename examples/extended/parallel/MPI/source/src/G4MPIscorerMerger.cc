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
#include "G4MPIscorerMerger.hh"
#include <map>
#include <ostream>
#include <algorithm>
#include <assert.h>
#include <functional>
#include "G4MPIutils.hh"

#define DMSG( LVL , MSG ) { if ( verbose > LVL ) { G4cout << MSG << G4endl; } }

namespace {
    //This class extends G4StatDouble to get access
    //to all data-members of the base class.
    //It provides two functions: pack and unpack
    //into a buffer for sending/receiving via MPI
    struct MPIStatDouble : public G4StatDouble {
      G4int verbose;
      inline void Pack(void* buffer,int bufferSize,int* position , MPI::Intracomm& comm ) const
      {
        DMSG(4,"Packing G4StatDouble(n,scale,sum_w,sum_w2,sum_wx,sum_wx2): "
             <<m_n<<" "<<m_scale<<" "<<m_sum_w<<" "<<m_sum_w2
             <<" "<<m_sum_wx<<" "<<m_sum_wx2);
        MPI_Pack(&m_n,1,MPI::INT,buffer,bufferSize,position,comm);
        const G4double data[]{m_scale,m_sum_w,m_sum_w2,m_sum_wx,m_sum_wx2};
        MPI_Pack(&data,5,MPI::DOUBLE,buffer,bufferSize,position,comm);
      }
      inline void UnPack(void* buffer,int bufferSize,int* position , MPI::Intracomm& comm ) {
        MPI_Unpack(buffer,bufferSize,position,&m_n,1,MPI::INT,comm);
        G4double data[5];
        MPI_Unpack(buffer,bufferSize,position,data,5,MPI::DOUBLE,comm);
        m_scale = data[0];
        m_sum_w = data[1];
        m_sum_w2= data[2];
        m_sum_wx= data[3];
        m_sum_wx2=data[4];
        DMSG(4,"UnPacking G4StatDouble(n,scale,sum_w,sum_w2,sum_wx,sum_wx2): "
             <<m_n<<" "<<m_scale<<" "<<m_sum_w<<" "<<m_sum_w2
             <<" "<<m_sum_wx<<" "<<m_sum_wx2);
      }
      MPIStatDouble(G4int ver = 0) : verbose(ver) {}
      MPIStatDouble(const G4StatDouble& rhs , G4int ver) : verbose(ver)
      {
        G4StatDouble::operator=(rhs);
      }
    };
}

G4MPIscorerMerger::G4MPIscorerMerger() :
  outputBuffer(nullptr),outputBufferSize(0),outputBufferPosition(0),bytesSent(0),
  ownsBuffer(false),scoringManager(nullptr),commSize(0),
  destinationRank(G4MPImanager::kRANK_MASTER),verbose(0)
{}

G4MPIscorerMerger::G4MPIscorerMerger(G4ScoringManager* mgr,
    G4int destination,
    G4int verbosity) :
    outputBuffer(nullptr),outputBufferSize(0),outputBufferPosition(0),bytesSent(0),
    ownsBuffer(false),
  scoringManager(mgr), commSize(0), destinationRank(destination),
      verbose(verbosity)
{
}

G4MPIscorerMerger::~G4MPIscorerMerger() {
  if ( ownsBuffer ) delete[] outputBuffer;
}

/* Format of the message.
 *
 * Input:
 * A vector of G4VScoringMesh, each of that is a
 * std::map<name:G4String,G4THitsMap<G4double>*> where
 * G4THitsMap<T> = std::map<int,T*>
 *
 * Output:
 * A buffer:
 *      [0] : numMesh : int             (**Begin Message**)
 *      [1] : meshID : int                (** Begin Mesh**)
 *      [2] : numMaps : int
 *      [3] : sizeName : int             (** Begin Map **)
 *      [4] : name[0] : char
 *      ...
 *      [...] : name[sizeName-1] : chare
 *      [...] : mapSize : int
 *      [...] : THitsMap.keys()[0] : int
 *      ...
 *      [...] : THitsMap.keys()[mapSize-1] : int
 *      [...] : THitsMap.values()[0] : double
 *      ...
 *      [...] : THitsMap.values()[mapSize-1] : double  (**End Map**)
 *      [...] : Next Map : repeat from (**Begin Map**)
 *      ...
 *      [...] : Next Mesh : repeat from (**Begin Mesh**)
 *
 *
 */

void G4MPIscorerMerger::Merge() {
  DMSG(0, "G4MPIscorerMerger::Merge called");
  const unsigned int myrank = G4MPImanager::GetManager()->GetRank();
  commSize = G4MPImanager::GetManager()->GetActiveSize();
  if ( commSize == 1 ) {
      DMSG(1,"Comm world size is 1, nothing to do");
      return;
  }
  const MPI::Intracomm* parentComm = G4MPImanager::GetManager()->GetComm();
  comm = parentComm->Dup();
  DestroyBuffer();

  //ANDREA:->
//  G4cout<<"Before sending: "<<G4endl;
//  scoringManager->GetMesh(0)->Dump();
//  for ( int i = 0 ; i < scoringManager->GetNumberOfMesh() ; ++i )  {
//  for ( auto e : scoringManager->GetMesh(i)->GetScoreMap() )
//    {
//      G4cout<<e.first<<" : "<<e.second<<G4endl;
//      for ( auto c: *(e.second->GetMap()) ) {
//          G4cout<<c.first<<"="<<*c.second<<G4endl;
//
//     }
//    }
//  }
  //ANDREA:<-

  bytesSent=0;
  const G4double sttime = MPI::Wtime();

  //Use G4MPIutils to optimize communications between ranks
  typedef std::function<void(unsigned int)> handler_t;
  using std::placeholders::_1;
  handler_t sender = std::bind(&G4MPIscorerMerger::Send , this , _1);
  handler_t receiver = std::bind(&G4MPIscorerMerger::Receive, this, _1);
  std::function<void(void)> barrier = std::bind(&MPI::Intracomm::Barrier,&comm);
  G4mpi::Merge( sender , receiver , barrier , commSize , myrank );

  //OLD Style p2p communications
  /*
  if ( myrank != destinationRank ) {
      DMSG(1,"Comm world size: "<<commSize<<" this rank is: "
          <<myrank<<" sending to rank "<<destinationRank
          <<" Number of mesh: "<< scoringManager->GetNumberOfMesh() );
      Send(destinationRank);
  } else {
      DMSG(1,"Comm world size: "<<commSize<<" this rank is: "
          <<myrank<<" receiving "
          <<" Number of mesh: "<< scoringManager->GetNumberOfMesh() );
      for ( unsigned int i = 0 ; i < commSize ; ++i ) {
          if ( i != myrank ) Receive(i);
      }
  }
*/
  const G4double elapsed = MPI::Wtime() - sttime;
  long total=0;
  comm.Reduce(&bytesSent,&total,1,MPI::LONG,MPI::SUM,destinationRank);
  if ( verbose > 0 && myrank == destinationRank ) {
      //Collect from ranks how much data was sent around
      G4cout<<"G4MPIscorerMerger::Merge() -data transfer performances: "
                  <<double(total)/1000./elapsed<<" kB/s"
                  <<" (Total Data Transfer= "<<double(total)/1000.<<" kB in "
                  <<elapsed<<" s)."<<G4endl;
  }
  //ANDREA:->
//  G4cout<<"After Receiving: "<<G4endl;
//  scoringManager->GetMesh(0)->Dump();
//  for ( int i = 0 ; i < scoringManager->GetNumberOfMesh() ; ++i )  {
//  for ( auto e : scoringManager->GetMesh(i)->GetScoreMap() )
//    {
//      G4cout<<e.first<<" : "<<e.second<<G4endl;
//      for ( auto c: *(e.second->GetMap()) ) {
//          G4cout<<c.first<<"="<<*c.second<<" (=2x"<<.5*(*c.second)<<")"<<G4endl;
//
//     }
//    }
//  }
  //ANDREA:<-
  comm.Free();
  DMSG(0,"G4MPIscorerMerger::Merge done.");
}

void G4MPIscorerMerger::Receive(const unsigned int source) {
  DMSG(1,"Receiving scorers");
 // DestroyBuffer();
      DMSG(2,"Receiving from: "<<source);
      MPI::Status status;
      comm.Probe(source, G4MPImanager::kTAG_CMDSCR, status);
      const G4int newbuffsize = status.Get_count(MPI::PACKED);
      DMSG(2,"Preparing to receive buffer of size: "<<newbuffsize);
      char* buffer = outputBuffer;
      if ( newbuffsize > outputBufferSize ) {
          DMSG(3,"New larger buffer expected, resize");
          //New larger buffer incoming, recreate buffer
          //TODO: use realloc?
          delete[] outputBuffer;
          buffer = new char[newbuffsize];
          //Avoid complains from valgrind (i'm not really sure why this is needed, but, beside the
          //small cpu penalty, we can live with that).)
          std::fill( buffer , buffer + newbuffsize , 0 );
          ownsBuffer = true;
      }
      SetupOutputBuffer(buffer,newbuffsize,0);
      comm.Recv(buffer, newbuffsize, MPI::PACKED, source,
          G4MPImanager::kTAG_CMDSCR, status);
      DMSG(3,"Buffer Size: "<<outputBufferSize<< " bytes at: "<<(void*)outputBuffer);
      UnPackAndMerge(scoringManager);
  DMSG(1,"Receiving of comamnd line scorers done");
}

void G4MPIscorerMerger::Send(const unsigned int destination) {
  DMSG(1,"Sending scorers "<<this);
  //Step 1: Setup buffer to pack/unpack data
  const G4int newbuffsize = CalculatePackSize(scoringManager);
  //DestroyBuffer();
  char* buffer = outputBuffer;
  if ( newbuffsize > outputBufferSize ) {
      delete[] outputBuffer;
      buffer = new char[newbuffsize];
      //Avoid complains from valgrind (i'm not really sure why this is needed, but, beside the
      //small cpu penalty, we can live with that).)
      std::fill( buffer , buffer+newbuffsize,0);
      ownsBuffer = true;
  }
  SetupOutputBuffer(buffer,newbuffsize,0);
  DMSG(3,"Buffer Size: "<<newbuffsize<< " bytes at: "<<(void*)outputBuffer);
  Pack(scoringManager);
  assert(outputBufferSize==outputBufferPosition);

  //Version 1: p2p communication
  comm.Send(outputBuffer, outputBufferSize, MPI::PACKED, destination, G4MPImanager::kTAG_CMDSCR);
  bytesSent += newbuffsize;
  //Receiver should use probe to get size of the package being sent
  DMSG(1,"Sending done");
}

void G4MPIscorerMerger::Pack(const G4ScoringManager* sm) {
  assert(sm!=nullptr);
  if ( outputBuffer == nullptr || outputBufferPosition>=outputBufferSize) {
      G4Exception("G4MPIscorerMerger::Pack(const G4ScoringManager*)",
                "MPI001",FatalException,
                "Call SetOututBuffer before trying to pack");
      return;
  }
  DMSG(2,"Starting packing of meshes, # meshes: "<<sm->GetNumberOfMesh());
  /*const*/ size_t numMeshes=sm->GetNumberOfMesh();//TODO: OLD MPI interface
  MPI_Pack(&numMeshes,1,MPI::UNSIGNED,
                outputBuffer,outputBufferSize,
                &outputBufferPosition,
                comm);
  for (size_t i = 0; i <numMeshes; ++i)
    {
      MPI_Pack(&i,1,MPI::UNSIGNED,
          outputBuffer,outputBufferSize,
          &outputBufferPosition,comm);
      Pack(sm->GetMesh(i));
    }
}

void G4MPIscorerMerger::UnPackAndMerge(const G4ScoringManager* sm) {
  assert(sm!=nullptr);
  if ( outputBuffer == nullptr || outputBufferPosition>=outputBufferSize) {
      G4Exception("G4MPIscorerMerger::UnPack(const G4ScroingManager*)",
                "MPI001",FatalException,
                "Call SetOututBuffer before trying to un-pack");
      return;
  }
  size_t numMeshes=0;
  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
      &numMeshes,1,MPI::UNSIGNED,comm);
  if ( numMeshes != sm->GetNumberOfMesh() ) {
      G4ExceptionDescription msg;
      msg << "Number of meshes to unpack ("<<numMeshes;
      msg <<") does not correspond to expected number ("<<sm->GetNumberOfMesh();
      msg<<")";
      G4Exception("G4MPIscorerMerger::UnPack(const G4ScroingManager*)",
                "MPI001",FatalException,msg);
      return;
  }

  size_t meshid=0;
  for ( size_t i = 0 ; i < numMeshes ; ++i ) {
      MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
          &meshid,1,MPI::UNSIGNED,comm);
      if ( meshid != i )  {
          G4ExceptionDescription msg;
          msg<<"Cannot unpack: expecting mesh "<<i<<" and found "<<meshid;
          msg<<" during unpack.";
          G4Exception("G4MPIscorerMerger::UnPack(const G4ScroingManager*)",
                    "MPI001",FatalException,msg);
          return;
      }
      G4VScoringMesh* original = sm->GetMesh(i);
      UnPackAndMerge(original);
  }
}

void G4MPIscorerMerger::Pack(const G4VScoringMesh* mesh) {
  assert(mesh!=nullptr);
  assert(outputBuffer!=nullptr);
  assert(outputBufferPosition<=outputBufferSize);
  DMSG(3,"Packing mesh: "<<mesh);

  auto map = mesh->GetScoreMap();
  /*const*/ size_t nummaps = map.size();//TODO: old MPI interface
  MPI_Pack(&nummaps,1,MPI::UNSIGNED,
      outputBuffer,outputBufferSize,
      &outputBufferPosition,comm);
  for ( const auto& ele: map ) {
      const G4String& name = ele.first;
      /*const*/ size_t ss = name.size();
      MPI_Pack(&ss,1,MPI::UNSIGNED,
          outputBuffer,outputBufferSize,
          &outputBufferPosition,comm);
#ifdef G4MPI_USE_MPI_PACK_NOT_CONST
      char* nn = new char[name.length()];
      std::copy(name.begin(),name.end(),nn); 
#else
      const char* nn = name.c_str();
#endif
      MPI_Pack(nn,ss,MPI::CHAR,outputBuffer,outputBufferSize,&outputBufferPosition,comm);
      Pack(ele.second);
#ifdef G4MPI_USE_MPI_PACK_NOT_CONST
      delete[] nn;
#endif
  }
}

void G4MPIscorerMerger::UnPackAndMerge(G4VScoringMesh* inmesh) {
  assert(outputBuffer!=nullptr);
  assert(outputBufferPosition<=outputBufferSize);
  assert(inmesh!=nullptr);
  DMSG(3,"Preparing to unpack a mesh and merge into: "<<inmesh);
  const G4String& detName = inmesh->GetWorldName();
  size_t nummaps = 0;
  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
      &nummaps,1,MPI::UNSIGNED,comm);
  for ( size_t i = 0 ; i < nummaps ; ++i ) {
      size_t nameSize = 0;
      MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
          &nameSize,1,MPI::UNSIGNED,comm);
      //Create a null-terminated c-string: needed later when converting this to a G4String
      //(Not sure: but issue reported by valgrind with the use of MPI_Unpack)
      char* name = new char[nameSize+1];
      std::fill(name,name+nameSize+1,0);
      MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
          name,nameSize,MPI::CHAR,comm);
      const G4String colname(name,nameSize);
      delete[] name;
      //This memory churn is very inefficient, but we cannot reuse the HitMap
      //because we cannot change the names
      //TODO: Evaluate change in HitMap class to allow for change of names
      //HitMap* hm = UnPackHitMap(detName,colname);
      HitStatDoubleMap* hm = UnPackHitStatDoubleMap(detName,colname);
      inmesh->Accumulate(hm);
      delete hm;
  }
}


//void G4MPIscorerMerger::Pack(const HitMap* sm) {
//  assert(sm!=nullptr);
//  assert(outputBuffer!=nullptr);
//  assert(outputBufferPosition<=outputBufferSize);
//  DMSG(3,"Packing hitmap: "<<sm<<" with: "<<sm->GetSize()<<" elements.");
//  /*const*/ size_t numEl = sm->GetSize();//TODO: old MPI implementation
//  MPI_Pack(&numEl,1,MPI::UNSIGNED,
//      outputBuffer,outputBufferSize,
//      &outputBufferPosition,comm);
//  const auto& theMap = *sm->GetMap();
//  std::vector<G4int> ids;
//  std::vector<G4double> vals;
//  std::transform(theMap.begin(),theMap.end(),std::back_inserter(ids),
//      [](decltype(*theMap.begin())& e){ return e.first;});
//  std::transform(theMap.begin(),theMap.end(),std::back_inserter(vals),
//      [](decltype(*theMap.begin())& e){ return *e.second;});
//  assert(ids.size()==vals.size()&&ids.size()==numEl);
//  MPI_Pack(ids.data(),ids.size(),MPI::INT,
//      outputBuffer,outputBufferSize,
//      &outputBufferPosition,comm);
//  MPI_Pack(vals.data(),vals.size(),MPI::DOUBLE,
//      outputBuffer,outputBufferSize,
//      &outputBufferPosition,comm);
//}

void G4MPIscorerMerger::Pack(const HitStatDoubleMap* sm) {
  assert(sm!=nullptr);
  assert(outputBuffer!=nullptr);
  assert(outputBufferPosition<=outputBufferSize);
  DMSG(3,"Packing hitmap: "<<sm<<" with: "<<sm->GetSize()<<" elements.");
  /*const*/ size_t numEl = sm->GetSize();//TODO: old MPI implementation
  MPI_Pack(&numEl,1,MPI::UNSIGNED,
      outputBuffer,outputBufferSize,
      &outputBufferPosition,comm);
  const auto& theMap = *sm->GetMap();
  std::vector<G4int> ids;
  std::transform(theMap.begin(),theMap.end(),std::back_inserter(ids),
      [](decltype(*theMap.begin())& e){ return e.first;});
  assert(/*ids.size()==vals.size()&&*/ids.size()==numEl);
  MPI_Pack(ids.data(),ids.size(),MPI::INT,outputBuffer,outputBufferSize,
      &outputBufferPosition,comm);
  for( const auto& e : theMap) {
     const MPIStatDouble sd(*e.second,verbose);
     sd.Pack(outputBuffer,outputBufferSize,&outputBufferPosition,comm);
  }
}

//HitMap* G4MPIscorerMerger::UnPackHitMap(const G4String& detName,
//                                        const G4String& colName) {
//  assert(outputBuffer!=nullptr);
//  assert(outputBufferPosition<=outputBufferSize);
//  DMSG(3,"Preparing to unpack a hit map for: "<<detName<<","<<colName);
//  size_t numEl =0 ;
//  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
//             &numEl,1,MPI::UNSIGNED,comm);
//  G4int* ids = new G4int[numEl];
//  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
//             ids,numEl,MPI::INT,comm);
//  G4double* vals = new G4double[numEl];
//  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
//      vals,numEl,MPI::DOUBLE,comm);
//  HitMap* result = new HitMap(detName,colName);
//  for ( unsigned int i = 0; i<numEl;++i) result->set(ids[i],vals[i]);
//  delete[] ids;
//  delete[] vals;
//  return result;
//}

HitStatDoubleMap* G4MPIscorerMerger::UnPackHitStatDoubleMap(
    const G4String& detName, const G4String& colName)
{
  assert(outputBuffer!=nullptr);
  assert(outputBufferPosition<=outputBufferSize);
  DMSG(3,"Preparing to unpack a hit map for: "<<detName<<","<<colName);
  size_t numEl =0 ;
  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
             &numEl,1,MPI::UNSIGNED,comm);
  DMSG(3,"Will receive "<<numEl<<" values");
  G4int* ids = new G4int[numEl];
  MPI_Unpack(outputBuffer,outputBufferSize,&outputBufferPosition,
             ids,numEl,MPI::INT,comm);
  HitStatDoubleMap* result = new HitStatDoubleMap(detName,colName);
  for ( unsigned int i = 0; i<numEl;++i) {
      MPIStatDouble sd(verbose);
      sd.UnPack(outputBuffer,outputBufferSize,&outputBufferPosition,comm);
      result->set(ids[i],sd);
    }
  delete[] ids;
  return result;
}


G4int G4MPIscorerMerger::CalculatePackSize(const G4ScoringManager* sm) const
{
  DMSG(3,"Calculating dimension of data to send");
  if ( sm == nullptr ) return 0;
  //Calcualte how much data each call to Pack* appends to the buffer
  //e.g. sizeof(data)
  //The number of sizeof here should match the number of calls to MPI_Pack

  //Pack(ScoringMgr)
  G4int size = sizeof(unsigned int);
  DMSG(3,"There are "<<sm->GetNumberOfMesh()<<" meshes.");
  //Loop on mesh
  for ( size_t i = 0 ; i<sm->GetNumberOfMesh() ; ++i ) {
      size += sizeof(unsigned int);//ID
      size += CalculatePackSize(sm->GetMesh(i));
  }
  return size;
}

G4int G4MPIscorerMerger::CalculatePackSize(const G4VScoringMesh* mesh) const
{
  DMSG(3,"Calculating size for mesh: "<<mesh);
  //PackSingleMesh(Mesh)
  G4int size = sizeof(unsigned int);//num maps
  auto map = mesh->GetScoreMap();
  for (const auto& ele : map ) {
      //PackHitsMap
      size += sizeof(unsigned int);//name size
      const G4String& name = ele.first;
      size += sizeof(char)*name.size();//name
      size += CalculatePackSize(ele.second);
  }
  DMSG(3,"mesh "<<mesh<<" size: "<<size);
  return size;
}

//G4int G4MPIscorerMerger::CalculatePackSize(const HitMap* map) const {
//  const G4int numEls = map->GetSize();
//  G4int size = sizeof(unsigned int);
//  size += sizeof(G4int)*numEls;
//  size += sizeof(G4double)*numEls;
//  DMSG(3,"HitMap "<<map<<" size: "<<size<<" in "<<numEls<<" elements.");
//  return size;
//}

G4int G4MPIscorerMerger::CalculatePackSize(const HitStatDoubleMap* map) const {
  const G4int numEls = map->GetSize();
  G4int size = sizeof(unsigned int);
  size += sizeof(G4int)*numEls;
  //G4StatDouble: 5 doubles and 1 int
  //Can I use sizeof(G4StatDouble)? NO sizeof(G4StatDouble)==56
  size += numEls*(sizeof(G4double)*5+sizeof(G4int));
  DMSG(3,"HitStatDoubleMap "<<map<<" size: "<<size<<" in "<<numEls<<" elements.");
  return size;
}

