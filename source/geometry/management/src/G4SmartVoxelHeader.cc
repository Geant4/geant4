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
// class G4SmartVoxelHeader implementation
//
// Define G4GEOMETRY_VOXELDEBUG for debugging information on G4cout
//
// 29.04.02 Use 3D voxelisation for non consuming replication - G.C.
// 18.04.01 Migrated to STL vector - G.C.
// 12.02.99 Introduction of new quality/smartless: max for (slices/cand) - S.G.
// 11.02.99 Voxels at lower levels are now built for collapsed slices - S.G.
// 21.07.95 Full implementation, supporting non divided physical volumes - P.K.
// 14.07.95 Initial version - stubb definitions only - P.K.
// --------------------------------------------------------------------

#include "G4SmartVoxelHeader.hh"

#include "G4ios.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"

#include "voxeldefs.hh"
#include "G4AffineTransform.hh"
#include "G4VSolid.hh"
#include "G4VPVParameterisation.hh"

// ***************************************************************************
// Constructor for topmost header, to begin voxel construction at a
// given logical volume.
// Constructs target List of volumes, calls "Build and refine" constructor.
// Assumes all daughters represent single volumes (ie. no divisions
// or parametric)
// ***************************************************************************
//
G4SmartVoxelHeader::G4SmartVoxelHeader(G4LogicalVolume* pVolume,
                                       G4int pSlice)
  : fminEquivalent(pSlice),
    fmaxEquivalent(pSlice),
    fparamAxis(kUndefined)
{
  std::size_t nDaughters = pVolume->GetNoDaughters();

  // Determine whether daughter is replicated
  //
  if ((nDaughters!=1) || (!pVolume->GetDaughter(0)->IsReplicated()))
  {
    // Daughter not replicated => conventional voxel Build
    // where each daughters extents are computed
    //
    BuildVoxels(pVolume);
  }
  else
  {
    // Single replicated daughter
    //
    BuildReplicaVoxels(pVolume);
  }
}

// ***************************************************************************
// Protected constructor:
// builds and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'. `pSlice' is used to set max
// and min equivalent slice nos for the header - they apply to the level
// of the header, not its nodes.
// ***************************************************************************
//
G4SmartVoxelHeader::G4SmartVoxelHeader(G4LogicalVolume* pVolume,
                                 const G4VoxelLimits& pLimits,
                                 const G4VolumeNosVector* pCandidates,
                                       G4int pSlice)
  : fminEquivalent(pSlice),
    fmaxEquivalent(pSlice),
    fparamAxis(kUndefined)
{
#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "**** G4SmartVoxelHeader::G4SmartVoxelHeader" << G4endl
         << "     Limits " << pLimits << G4endl
         << "     Candidate #s = " ;
  for (auto i=0; i<pCandidates->size(); ++i)
  {
    G4cout << (*pCandidates)[i] << " ";
  }
  G4cout << G4endl;
#endif   

  BuildVoxelsWithinLimits(pVolume,pLimits,pCandidates);
}

// ***************************************************************************
// Destructor:
// deletes all proxies and underlying objects.
// ***************************************************************************
//
G4SmartVoxelHeader::~G4SmartVoxelHeader()
{
  // Manually destroy underlying nodes/headers
  // Delete collected headers and nodes once only
  //
  std::size_t node, proxy, maxNode=fslices.size();
  G4SmartVoxelProxy* lastProxy = nullptr;
  G4SmartVoxelNode *dyingNode, *lastNode=nullptr;
  G4SmartVoxelHeader *dyingHeader, *lastHeader=nullptr;

  for (node=0; node<maxNode; ++node)
  {
    if (fslices[node]->IsHeader())
    {
      dyingHeader = fslices[node]->GetHeader();
      if (lastHeader != dyingHeader)
      {
        lastHeader = dyingHeader;
        lastNode = nullptr;
        delete dyingHeader;
      }
    }
    else
    {
      dyingNode = fslices[node]->GetNode();
      if (dyingNode != lastNode)
      {
        lastNode = dyingNode;
        lastHeader = nullptr;
        delete dyingNode;
      }
    }
  }
  // Delete proxies
  //
  for (proxy=0; proxy<maxNode; ++proxy)
  {
    if (fslices[proxy] != lastProxy)
    {
      lastProxy = fslices[proxy];
      delete lastProxy;
    }
  }
  // Don't need to clear slices
  // fslices.clear();
}

// ***************************************************************************
// Equality operator: returns true if contents are equivalent.
// Implies a deep search through contained nodes/header.
// Compares headers' axes,sizes,extents. Returns false if different.
// For each contained proxy, determines whether node/header, compares and
// returns if different. Compares and returns if proxied nodes/headers
// are different.
// ***************************************************************************
//
G4bool G4SmartVoxelHeader::operator == (const G4SmartVoxelHeader& pHead) const
{
  if ( (GetAxis()      == pHead.GetAxis())
    && (GetNoSlices()  == pHead.GetNoSlices())
    && (GetMinExtent() == pHead.GetMinExtent())
    && (GetMaxExtent() == pHead.GetMaxExtent()) )
  {
    std::size_t node, maxNode;
    G4SmartVoxelProxy *leftProxy, *rightProxy;
    G4SmartVoxelHeader *leftHeader, *rightHeader;
    G4SmartVoxelNode *leftNode, *rightNode;

    maxNode = GetNoSlices();
    for (node=0; node<maxNode; ++node)
    {
      leftProxy  = GetSlice(node);
      rightProxy = pHead.GetSlice(node);
      if (leftProxy->IsHeader())
      {
        if (rightProxy->IsNode())
        {
          return false;
        }
        else
        {
          leftHeader  = leftProxy->GetHeader();
          rightHeader = rightProxy->GetHeader();
          if (!(*leftHeader == *rightHeader))
          {
            return false;
          }
        }
      }
      else
      {
        if (rightProxy->IsHeader())
        {
          return false;
        }
        else
        {
          leftNode  = leftProxy->GetNode();
          rightNode = rightProxy->GetNode();
          if (!(*leftNode == *rightNode))
          {
            return false;
          }
        }
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

// ***************************************************************************
// Builds voxels for daughters specified volume, in NON-REPLICATED case
// o Create List of target volume nos (all daughters; 0->noDaughters-1)
// o BuildWithinLimits does Build & also determines mother dimensions.
// ***************************************************************************
//
void G4SmartVoxelHeader::BuildVoxels(G4LogicalVolume* pVolume)
{
  G4VoxelLimits limits;   // Create `unlimited' limits object
  std::size_t nDaughters = pVolume->GetNoDaughters();

  G4VolumeNosVector targetList;
  targetList.reserve(nDaughters);
  for (std::size_t i=0; i<nDaughters; ++i)
  {
    targetList.push_back((G4int)i);
  }
  BuildVoxelsWithinLimits(pVolume, limits, &targetList);
}

// ***************************************************************************
// Builds voxels for specified volume containing a single replicated volume.
// If axis is not specified (i.e. "kUndefined"), 3D voxelisation is applied,
// and the best axis is determined according to heuristics as for placements.
// ***************************************************************************
//
void G4SmartVoxelHeader::BuildReplicaVoxels(G4LogicalVolume* pVolume)
{
  G4VPhysicalVolume* pDaughter = nullptr;

  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;

  // Consistency check: pVolume should contain single replicated volume
  //
  if ( (pVolume->GetNoDaughters()==1)
    && (pVolume->GetDaughter(0)->IsReplicated()==true) )
  {
    // Obtain replication data
    //
    pDaughter = pVolume->GetDaughter(0);
    pDaughter->GetReplicationData(axis,nReplicas,width,offset,consuming);
    fparamAxis = axis;
    if ( consuming == false )
    {
      G4VoxelLimits limits;   // Create `unlimited' limits object
      G4VolumeNosVector targetList;
      targetList.reserve(nReplicas);
      for (auto i=0; i<nReplicas; ++i)
      {
        targetList.push_back(i);
      }
      if (axis != kUndefined)
      {
        // Apply voxelisation along the specified axis only

        G4ProxyVector* pSlices=BuildNodes(pVolume,limits,&targetList,axis);
        faxis = axis;
        fslices = *pSlices;
        delete pSlices;

        // Calculate and set min and max extents given our axis
        //
        const G4AffineTransform origin;
        pVolume->GetSolid()->CalculateExtent(faxis, limits, origin,
                                             fminExtent, fmaxExtent);
        // Calculate equivalent nos
        //
        BuildEquivalentSliceNos();
        CollectEquivalentNodes();   // Collect common nodes
      }
      else
      {
        // Build voxels similarly as for normal placements considering
        // all three cartesian axes.

        BuildVoxelsWithinLimits(pVolume, limits, &targetList);
      }
    }
    else
    {
      // Replication is consuming -> Build voxels directly
      //
      // o Cartesian axes - range is -width*nREplicas/2 to +width*nREplicas/2
      //                    nReplicas replications result
      // o Radial axis (rho) = range is 0 to width*nReplicas
      //                    nReplicas replications result
      // o Phi axi       - range is offset to offset+width*nReplicas radians
      //
      // Equivalent slices no computation & collection not required - all
      // slices are different
      //
      switch (axis)
      {
        case kXAxis:
        case kYAxis:
        case kZAxis:
          fminExtent = -width*nReplicas*0.5;
          fmaxExtent =  width*nReplicas*0.5;
          break;
        case kRho:
          fminExtent = offset;
          fmaxExtent = width*nReplicas+offset;
          break;
        case kPhi:
          fminExtent = offset;
          fmaxExtent = offset+width*nReplicas;
          break;
        default:
          G4Exception("G4SmartVoxelHeader::BuildReplicaVoxels()",
                      "GeomMgt0002", FatalException, "Illegal axis.");
          break;
      }  
      faxis = axis;   // Set axis
      BuildConsumedNodes(nReplicas);
      if ( (axis==kXAxis) || (axis==kYAxis) || (axis==kZAxis) )
      {
        // Sanity check on extent
        //
        G4double emin = kInfinity, emax = -kInfinity;
        G4VoxelLimits limits;
        G4AffineTransform origin;
        pVolume->GetSolid()->CalculateExtent(axis, limits, origin, emin, emax);
        if ( (std::fabs((emin-fminExtent)/fminExtent) +
              std::fabs((emax-fmaxExtent)/fmaxExtent)) > 0.05)
        {
          std::ostringstream message;
          message << "Sanity check: wrong solid extent." << G4endl
                  << "        Replicated geometry, logical volume: "
                  << pVolume->GetName();
          G4Exception("G4SmartVoxelHeader::BuildReplicaVoxels",
                      "GeomMgt0002", FatalException, message);
        }
      }
    }
  }
  else
  {
    G4Exception("G4SmartVoxelHeader::BuildReplicaVoxels", "GeomMgt0002",
                FatalException, "Only one replicated daughter is allowed !");
  }
}

// ***************************************************************************
// Builds `consumed nodes': nReplicas nodes each containing one replication,
// numbered in sequence 0->nReplicas-1
// o Modifies fslices `in place'
// o faxis,fminExtent,fmaxExtent NOT modified.
// ***************************************************************************
//
void G4SmartVoxelHeader::BuildConsumedNodes(G4int nReplicas)
{
  G4int nNode, nVol;
  G4SmartVoxelNode* pNode;
  G4SmartVoxelProxy* pProxyNode;

  // Create and fill nodes in temporary G4NodeVector (on stack)
  //
  G4NodeVector nodeList;
  nodeList.reserve(nReplicas);
  for (nNode=0; nNode<nReplicas; ++nNode)
  {
    pNode = new G4SmartVoxelNode(nNode);
    if (pNode == nullptr)
    {
      G4Exception("G4SmartVoxelHeader::BuildConsumedNodes()", "GeomMgt0003",
                  FatalException, "Node allocation error.");
    }
    nodeList.push_back(pNode);
  }
  for (nVol=0; nVol<nReplicas; ++nVol)
  {
    nodeList[nVol]->Insert(nVol);   // Insert replication of number
  }                                 // identical to voxel number

  // Create & fill proxy List `in place' by modifying instance data fslices
  //
  fslices.clear();
  for (nNode=0; nNode<nReplicas; ++nNode)
  {
    pProxyNode = new G4SmartVoxelProxy(nodeList[nNode]);
    if (!pProxyNode)
    {
      G4Exception("G4SmartVoxelHeader::BuildConsumedNodes()", "GeomMgt0003",
                  FatalException, "Proxy node allocation error.");
    }
    fslices.push_back(pProxyNode);
  }
}

// ***************************************************************************
// Builds and refines voxels between specified limits, considering only
// the physical volumes numbered `pCandidates'.
// o Chooses axis
// o Determines min and max extents (of mother solid) within limits.
// ***************************************************************************
//
void
G4SmartVoxelHeader::BuildVoxelsWithinLimits(G4LogicalVolume* pVolume,
                                            G4VoxelLimits pLimits,
                                      const G4VolumeNosVector* pCandidates)
{
  // Choose best axis for slicing by:
  // 1. Trying all unlimited cartesian axes
  // 2. Select axis which gives greatest no slices

  G4ProxyVector *pGoodSlices=nullptr, *pTestSlices, *tmpSlices;
  G4double goodSliceScore=kInfinity, testSliceScore;
  EAxis goodSliceAxis = kXAxis;
  EAxis testAxis      = kXAxis;
  std::size_t node, maxNode, iaxis;
  G4VoxelLimits noLimits;

  // Try all non-limited cartesian axes
  //
  for (iaxis=0; iaxis<3; ++iaxis)
  {
    switch(iaxis)
    {
      case 0:
        testAxis = kXAxis;
        break;
      case 1:
        testAxis = kYAxis;
        break;
      case 2:
        testAxis = kZAxis;
        break;
    }
    if (!pLimits.IsLimited(testAxis))
    {
      pTestSlices = BuildNodes(pVolume,pLimits,pCandidates,testAxis);
      testSliceScore = CalculateQuality(pTestSlices);
      if ( (!pGoodSlices) || (testSliceScore<goodSliceScore) )
      {
        goodSliceAxis  = testAxis;
        goodSliceScore = testSliceScore;
        tmpSlices      = pGoodSlices;
        pGoodSlices    = pTestSlices;
        pTestSlices    = tmpSlices;
      }
      if (pTestSlices)
      {
        // Destroy pTestSlices and all its contents
        //
        maxNode=pTestSlices->size();
        for (node=0; node<maxNode; ++node)
        {
          delete (*pTestSlices)[node]->GetNode();
        }
        G4SmartVoxelProxy* tmpProx;
        while (pTestSlices->size()>0)  // Loop checking, 06.08.2015, G.Cosmo
        {
          tmpProx = pTestSlices->back();
          pTestSlices->pop_back();
          for (auto i=pTestSlices->cbegin(); i!=pTestSlices->cend(); )
          {
            if (*i==tmpProx)
            {
              i = pTestSlices->erase(i);
            }
            else
            {
              ++i;
            }
          }
          if ( tmpProx ) { delete tmpProx; }
        }
        delete pTestSlices;
      }
    }
  }
  // Check for error case.. when limits already 3d,
  // so cannot select a new axis
  //
  if (!pGoodSlices)
  {
    G4Exception("G4SmartVoxelHeader::BuildVoxelsWithinLimits()",
                "GeomMgt0002", FatalException,
                "Cannot select more than 3 axis for optimisation.");
    return;
  }

  // 
  // We have selected pGoodSlices, with a score testSliceScore
  //

  // Store chosen axis, slice ptr
  //
  fslices =* pGoodSlices; // Set slice information, copy ptrs in collection
  delete pGoodSlices;     // Destroy slices vector, but not contained
                          // proxies or nodes
  faxis = goodSliceAxis;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << G4endl << "     Volume = " << pVolume->GetName()
         << G4endl << "     Selected axis = " << faxis << G4endl;
  for (auto islice=0; islice<fslices.size(); ++islice)
  {
    G4cout << "     Node #" << islice << " = {";
    for (auto j=0; j<fslices[islice]->GetNode()->GetNoContained(); ++j)
    {
      G4cout << " " << fslices[islice]->GetNode()->GetVolume(j);
    }
    G4cout << " }" << G4endl;
  }
  G4cout << G4endl;
#endif

  // Calculate and set min and max extents given our axis
  //
  G4VSolid* outerSolid = pVolume->GetSolid();
  const G4AffineTransform origin;
  if(!outerSolid->CalculateExtent(faxis,pLimits,origin,fminExtent,fmaxExtent))
  {
    outerSolid->CalculateExtent(faxis,noLimits,origin,fminExtent,fmaxExtent);
  }

  // Calculate equivalent nos
  //
  BuildEquivalentSliceNos();
  CollectEquivalentNodes();      // Collect common nodes
  RefineNodes(pVolume, pLimits); // Refine nodes creating headers

  // No common headers can exist because collapsed by construction
}

// ***************************************************************************
// Calculates and stores the minimum and maximum equivalent neighbour
// values for all slices at our level.
//
// Precondition: all slices are nodes.
// For each potential start of a group of equivalent nodes:
// o searches forwards in fslices to find group end
// o loops from start to end setting start and end slices.
// ***************************************************************************
//
void G4SmartVoxelHeader::BuildEquivalentSliceNos()
{
  std::size_t sliceNo, minNo, maxNo, equivNo;
  std::size_t maxNode = fslices.size();
  G4SmartVoxelNode *startNode, *sampleNode;
  for (sliceNo=0; sliceNo<maxNode; ++sliceNo)
  {
    minNo = sliceNo;

    // Get first node (see preconditions - will throw exception if a header)
    //
    startNode = fslices[minNo]->GetNode();

    // Find max equivalent
    //
    for (equivNo=minNo+1; equivNo<maxNode; ++equivNo)
    {
      sampleNode = fslices[equivNo]->GetNode();
      if (!((*startNode) == (*sampleNode))) { break; }
    }
    maxNo = equivNo-1;
    if (maxNo != minNo)
    {
      // Set min and max nos
      //
      for (equivNo=minNo; equivNo<=maxNo; ++equivNo)
      {
        sampleNode = fslices[equivNo]->GetNode();
        sampleNode->SetMinEquivalentSliceNo((G4int)minNo);
        sampleNode->SetMaxEquivalentSliceNo((G4int)maxNo);
      }
      // Advance outer loop to end of equivalent group
      //
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Collects common nodes at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
//
// Preconditions:
// o the slices have not previously be "collected"
// o all of the slices are nodes.
// ***************************************************************************
//
void G4SmartVoxelHeader::CollectEquivalentNodes()
{
  std::size_t sliceNo, maxNo, equivNo;
  std::size_t maxNode=fslices.size();
  G4SmartVoxelNode* equivNode;
  G4SmartVoxelProxy* equivProxy;

  for (sliceNo=0; sliceNo<maxNode; ++sliceNo)
  {
    equivProxy=fslices[sliceNo];

    // Assumption (see preconditions): all slices are nodes
    //
    equivNode = equivProxy->GetNode();
    maxNo = equivNode->GetMaxEquivalentSliceNo();
    if (maxNo != sliceNo)
    {
#ifdef G4GEOMETRY_VOXELDEBUG
      G4cout << "**** G4SmartVoxelHeader::CollectEquivalentNodes" << G4endl
             << "     Collecting Nodes = " 
             << sliceNo << " - " << maxNo << G4endl;
#endif
      // Do collection between sliceNo and maxNo inclusive
      //
      for (equivNo=sliceNo+1; equivNo<=maxNo; ++equivNo)
      {
        delete fslices[equivNo]->GetNode();
        delete fslices[equivNo];
        fslices[equivNo] = equivProxy;
      }
      sliceNo = maxNo;
    }
  }
}

// ***************************************************************************
// Collects common headers at our level, deleting all but one to save
// memory, and adjusting stored slice pointers appropriately.
// 
// Preconditions:
// o if a header forms part of a range of equivalent slices
//   (ie. GetMaxEquivalentSliceNo()>GetMinEquivalentSliceNo()),
//   it is assumed that all slices in the range are headers.
// o this will be true if a constant Expression is used to evaluate
//   when to refine nodes.
// ***************************************************************************
//
void G4SmartVoxelHeader::CollectEquivalentHeaders()
{
  std::size_t sliceNo, maxNo, equivNo;
  std::size_t maxNode = fslices.size();
  G4SmartVoxelHeader *equivHeader, *sampleHeader;
  G4SmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; ++sliceNo)
  {
    equivProxy = fslices[sliceNo];
    if (equivProxy->IsHeader())
    {
      equivHeader = equivProxy->GetHeader();
      maxNo = equivHeader->GetMaxEquivalentSliceNo();
      if (maxNo != sliceNo)
      {
        // Attempt collection between sliceNo and maxNo inclusive:
        // look for common headers. All slices between sliceNo and maxNo
        // are guaranteed to be headers but may not have equal contents
        //
#ifdef G4GEOMETRY_VOXELDEBUG
        G4cout << "**** G4SmartVoxelHeader::CollectEquivalentHeaders" << G4endl
               << "     Collecting Headers =";
#endif
        for (equivNo=sliceNo+1; equivNo<=maxNo; ++equivNo)
        {
          sampleHeader = fslices[equivNo]->GetHeader();
          if ( (*sampleHeader) == (*equivHeader) )
          {
#ifdef G4GEOMETRY_VOXELDEBUG
            G4cout << " " << equivNo;
#endif
            // Delete sampleHeader + proxy and replace with equivHeader/Proxy
            //
            delete sampleHeader;
            delete fslices[equivNo];
            fslices[equivNo] = equivProxy;
          }
          else
          {
            // Not equal. Set this header to be
            // the current header for comparisons
            //
            equivProxy  = fslices[equivNo];
            equivHeader = equivProxy->GetHeader();
          }

        }
#ifdef G4GEOMETRY_VOXELDEBUG
        G4cout << G4endl;
#endif
        // Skip past examined slices
        //
        sliceNo = maxNo;
      }
    }
  }
}

// ***************************************************************************
// Builds the nodes corresponding to slices between the specified limits
// and along the specified axis, using candidate volume no.s in the vector
// pCandidates. If the `daughters' are replicated volumes (ie. the logical
// volume has a single replicated/parameterised volume for a daughter)
// the candidate no.s are interpreted as PARAMETERISED volume no.s & 
// PARAMETERISATIONs are applied to compute transformations & solid
// dimensions appropriately. The volume must be parameterised - ie. has a
// parameterisation object & non-consuming) - in this case.
// 
// Returns pointer to built node "structure" (guaranteed non NULL) consisting
// of G4SmartVoxelNodeProxies referring to G4SmartVoxelNodes.
// ***************************************************************************
//
G4ProxyVector* G4SmartVoxelHeader::BuildNodes(G4LogicalVolume* pVolume,
                                              G4VoxelLimits pLimits,
                                        const G4VolumeNosVector* pCandidates,
                                              EAxis pAxis)
{
  G4double motherMinExtent= kInfinity, motherMaxExtent= -kInfinity,
           targetMinExtent= kInfinity, targetMaxExtent= -kInfinity;
  G4VPhysicalVolume* pDaughter = nullptr;
  G4VPVParameterisation* pParam = nullptr;
  G4VSolid *targetSolid;
  G4AffineTransform targetTransform;
  G4bool replicated;
  std::size_t nCandidates = pCandidates->size();
  std::size_t nVol, nNode, targetVolNo;
  G4VoxelLimits noLimits;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "**** G4SmartVoxelHeader::BuildNodes" << G4endl
         << "     Limits = " << pLimits << G4endl
         << "       Axis = " << pAxis << G4endl
         << " Candidates = " << nCandidates << G4endl;
#endif

  // Compute extent of logical volume's solid along this axis
  // NOTE: results stored locally and not preserved/reused
  //
  G4VSolid* outerSolid = pVolume->GetSolid();
  const G4AffineTransform origin;
  if( !outerSolid->CalculateExtent(pAxis, pLimits, origin,
                                   motherMinExtent, motherMaxExtent) )
  {
    outerSolid->CalculateExtent(pAxis, noLimits, origin,
                                motherMinExtent, motherMaxExtent);
  }
  G4VolumeExtentVector minExtents(nCandidates,0.);
  G4VolumeExtentVector maxExtents(nCandidates,0.);

  if ( (pVolume->GetNoDaughters() == 1)
    && (pVolume->GetDaughter(0)->IsReplicated() == true) )
  {
    // Replication data not required: only parameterisation object 
    // and volume no. List used
    //
    pDaughter = pVolume->GetDaughter(0);
    pParam = pDaughter->GetParameterisation();
    if (pParam == nullptr)
    {
      std::ostringstream message;
      message << "PANIC! - Missing parameterisation." << G4endl
              << "         Replicated volume with no parameterisation object !";
      G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0003",
                  FatalException, message);
      return nullptr;
    }

    // Setup daughter's transformations
    //
    targetTransform = G4AffineTransform(pDaughter->GetRotation(),
                                        pDaughter->GetTranslation());
    replicated = true;
  }
    else
  {
    replicated = false;
  }
    
  // Compute extents
  //
  for (nVol=0; nVol<nCandidates; ++nVol)
  {
    targetVolNo = (*pCandidates)[nVol];
    if (replicated == false)
    {
      pDaughter = pVolume->GetDaughter(targetVolNo);

      // Setup daughter's transformations
      //
      targetTransform = G4AffineTransform(pDaughter->GetRotation(),
                                          pDaughter->GetTranslation());
      // Get underlying (and setup) solid
      //
      targetSolid = pDaughter->GetLogicalVolume()->GetSolid();
    }
    else
    {
      // Find  solid
      //
      targetSolid = pParam->ComputeSolid((G4int)targetVolNo,pDaughter);

      // Setup solid
      //
      targetSolid->ComputeDimensions(pParam,(G4int)targetVolNo,pDaughter);

      // Setup transform
      //
      pParam->ComputeTransformation((G4int)targetVolNo,pDaughter);
      targetTransform = G4AffineTransform(pDaughter->GetRotation(),
                                          pDaughter->GetTranslation());
    }
    // Calculate extents
    //
    if(!targetSolid->CalculateExtent(pAxis, pLimits, targetTransform,
                                     targetMinExtent, targetMaxExtent))
    {
      targetSolid->CalculateExtent(pAxis, noLimits, targetTransform,
                                   targetMinExtent,targetMaxExtent);
    }
    minExtents[nVol] = targetMinExtent;
    maxExtents[nVol] = targetMaxExtent;

#ifdef G4GEOMETRY_VOXELDEBUG
   G4cout << "---------------------------------------------------" << G4endl
          << "     Volume = " << pDaughter->GetName() << G4endl
          << " Min Extent = " << targetMinExtent << G4endl
          << " Max Extent = " << targetMaxExtent << G4endl
          << "---------------------------------------------------" << G4endl;
#endif

    // Check not entirely outside mother when processing toplevel nodes
    //
    if ( (!pLimits.IsLimited()) && ((targetMaxExtent<=motherMinExtent)
                                  ||(targetMinExtent>=motherMaxExtent)) )
    {
      std::ostringstream message;
      message << "PANIC! - Overlapping daughter with mother volume." << G4endl
              << "         Daughter physical volume "
              << pDaughter->GetName() << G4endl
              << "         is entirely outside mother logical volume "
              << pVolume->GetName() << " !!";
      G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0002",
                  FatalException, message);
    }

#ifdef G4GEOMETRY_VOXELDEBUG
    // Check for straddling volumes when debugging.
    // If a volume is >kStraddlePercent percent over the mother
    // boundary, print a warning.
    //
    if (!pLimits.IsLimited())
    {
      G4double width;
      G4int kStraddlePercent = 5;
      width = maxExtents[nVol]-minExtents[nVol];
      if ( (((motherMinExtent-minExtents[nVol])*100/width) > kStraddlePercent)
         ||(((maxExtents[nVol]-motherMaxExtent)*100/width) > kStraddlePercent) )
      {
        G4cout << "**** G4SmartVoxelHeader::BuildNodes" << G4endl
               << "     WARNING : Daughter # " << nVol
               << " name = " << pDaughter->GetName() << G4endl
               << "     Crosses mother boundary of logical volume, name = " 
               << pVolume->GetName() << G4endl
               << "     by more than " << kStraddlePercent 
               << "%" << G4endl;
      }
    }
#endif
  }

  // Extents of all daughters known

  // Calculate minimum slice width, only including volumes inside the limits
  //
  G4double minWidth = kInfinity;
  G4double currentWidth;
  for (nVol=0; nVol<nCandidates; ++nVol)
  {
    // currentWidth should -always- be a positive value. Inaccurate computed extent
    // from the solid or situations of malformed geometries (overlaps) may lead to
    // negative values and therefore unpredictable crashes !
    //
    currentWidth = std::abs(maxExtents[nVol]-minExtents[nVol]);
    if ( (currentWidth<minWidth)
      && (maxExtents[nVol]>=pLimits.GetMinExtent(pAxis))
      && (minExtents[nVol]<=pLimits.GetMaxExtent(pAxis)) )
    {
      minWidth = currentWidth;
    }
  }

  // No. of Nodes formula - nearest integer to
  // mother width/half min daughter width +1
  //
  G4double noNodesExactD = ((motherMaxExtent-motherMinExtent)*2.0/minWidth)+1.0;

  // Compare with "smartless quality", i.e. the average number of slices
  // used per contained volume.
  //
  G4double smartlessComputed = noNodesExactD / nCandidates;
  G4double smartlessUser = pVolume->GetSmartless();
  G4double smartless = (smartlessComputed <= smartlessUser)
                       ? smartlessComputed : smartlessUser;
  G4double noNodesSmart = smartless*nCandidates;
  G4int    noNodesExactI = G4int(noNodesSmart);
  G4long   noNodes = ((noNodesSmart-noNodesExactI)>=0.5)
                     ? noNodesExactI+1 : noNodesExactI;
  if( noNodes == 0 ) { noNodes=1; }

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "     Smartless computed = " << smartlessComputed << G4endl
         << "     Smartless volume = " << smartlessUser
         << " => # Smartless = " << smartless << G4endl;
  G4cout << "     Min width = " << minWidth
         << " => # Nodes = " << noNodes << G4endl;
#endif

  if (noNodes>kMaxVoxelNodes)
  {
    noNodes=kMaxVoxelNodes;
#ifdef G4GEOMETRY_VOXELDEBUG
    G4cout << "     Nodes Clipped to = " << kMaxVoxelNodes << G4endl;
#endif   
  }
  G4double nodeWidth = (motherMaxExtent-motherMinExtent)/noNodes;

  // Create G4VoxelNodes. Will Add proxies before setting fslices
  //
  auto* nodeList = new G4NodeVector();
  if (nodeList == nullptr)
  {
    G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0003",
                FatalException, "NodeList allocation error.");
    return nullptr;
  }
  nodeList->reserve(noNodes);

  for (nNode=0; G4long(nNode)<noNodes; ++nNode)
  {
    G4SmartVoxelNode *pNode;
    pNode = new G4SmartVoxelNode((G4int)nNode);
    if (pNode == nullptr)
    {
      G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0003",
                  FatalException, "Node allocation error.");
      return nullptr;
    }
    nodeList->push_back(pNode);
  }

  // All nodes created (empty)

  // Fill nodes: Step through extent lists
  //
  for (nVol=0; nVol<nCandidates; ++nVol)
  {
    G4long nodeNo, minContainingNode, maxContainingNode;
    minContainingNode = (minExtents[nVol]-motherMinExtent)/nodeWidth;
    maxContainingNode = (maxExtents[nVol]-motherMinExtent)/nodeWidth;

    // Only add nodes that are inside the limits of the axis
    //
    if ( (maxContainingNode>=0) && (minContainingNode<noNodes) )
    {
      // If max extent is on max boundary => maxContainingNode=noNodes:
      // should be one less as nodeList has noNodes entries
      //
      if (maxContainingNode>=noNodes)
      {
        maxContainingNode = noNodes-1;
      }
      //
      // Protection against protruding volumes
      //
      if (minContainingNode<0)
      {
        minContainingNode = 0;
      }
      for (nodeNo=minContainingNode; nodeNo<=maxContainingNode; ++nodeNo)
      {
        (*nodeList)[nodeNo]->Insert((*pCandidates)[nVol]);
      }
    }
  }

  // All nodes filled

  // Create proxy List : caller has deletion responsibility
  // (but we must delete nodeList *itself* - not the contents)
  //
  auto* proxyList = new G4ProxyVector();
  if (proxyList == nullptr)
  {
    G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0003",
                FatalException, "Proxy list allocation error.");
    return nullptr;
  }
  proxyList->reserve(noNodes);

  //
  // Fill proxy List
  //
  for (nNode=0; G4long(nNode)<noNodes; ++nNode)
  {
    // Get rid of possible excess capacity in the internal node vector
    //
    ((*nodeList)[nNode])->Shrink();
    auto* pProxyNode = new G4SmartVoxelProxy((*nodeList)[nNode]);
    if (pProxyNode == nullptr)
    {
      G4Exception("G4SmartVoxelHeader::BuildNodes()", "GeomMgt0003",
                  FatalException, "Proxy node allocation failed.");
      return nullptr;
    }
    proxyList->push_back(pProxyNode);
  }
  delete nodeList;
  return proxyList;
}

// ***************************************************************************
// Calculate a "quality value" for the specified vector of voxels.
// The value returned should be >0 and such that the smaller the number
// the higher the quality of the slice.
//
// Preconditions: pSlice must consist of G4SmartVoxelNodeProxies only
// Process:
// o Examine each node in turn, summing:
//      no. of non-empty nodes
//      no. of volumes in each node
// o Calculate Quality=sigma(volumes in nod)/(no. of non-empty nodes)
//      if all nodes empty, return kInfinity
// o Call G4Exception on finding a G4SmartVoxelHeaderProxy
// ***************************************************************************
//
G4double G4SmartVoxelHeader::CalculateQuality(G4ProxyVector *pSlice)
{
  G4double quality;
  std::size_t nNodes = pSlice->size();
  std::size_t noContained, maxContained=0, sumContained=0, sumNonEmptyNodes=0;
  G4SmartVoxelNode *node;

  for (std::size_t i=0; i<nNodes; ++i)
  {
    if ((*pSlice)[i]->IsNode())
    {
      // Definitely a node. Add info to running totals
      //
      node = (*pSlice)[i]->GetNode();
      noContained = node->GetNoContained();
      if (noContained)
      {
        ++sumNonEmptyNodes;
        sumContained += noContained;
        //
        // Calc maxContained for statistics
        //
        if (noContained>maxContained)
        {
          maxContained = noContained;
        }
      }
    }
    else
    {
      G4Exception("G4SmartVoxelHeader::CalculateQuality()", "GeomMgt0001",
                  FatalException, "Not applicable to replicated volumes.");
    }
  }

  // Calculate quality with protection against no non-empty nodes
  //
  if (sumNonEmptyNodes)
  {
    quality = sumContained/sumNonEmptyNodes;
  }
  else
  {
    quality = kInfinity;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "**** G4SmartVoxelHeader::CalculateQuality" << G4endl
         << "     Quality = " << quality << G4endl
         << "     Nodes = " << nNodes 
         << " of which " << sumNonEmptyNodes << " non empty" << G4endl
         << "     Max Contained = " << maxContained << G4endl;
#endif

  return quality;
}

// ***************************************************************************
// Examined each contained node, refines (creates a replacement additional
// dimension of voxels) when there is more than one voxel in the slice.
// Does not refine further if already limited in two dimensions (=> this
// is the third level of limits)
//
// Preconditions: slices (nodes) have been built.
// ***************************************************************************
//
void G4SmartVoxelHeader::RefineNodes(G4LogicalVolume* pVolume,
                                     G4VoxelLimits pLimits)
{
  std::size_t refinedDepth=0, minVolumes;
  std::size_t maxNode = fslices.size();

  if (pLimits.IsXLimited()) 
  {
    ++refinedDepth;
  }
  if (pLimits.IsYLimited()) 
  {
    ++refinedDepth;
  }
  if (pLimits.IsZLimited()) 
  {
    ++refinedDepth;
  }

  // Calculate minimum number of volumes necessary to refine
  //
  switch (refinedDepth)
  {
    case 0:
      minVolumes=kMinVoxelVolumesLevel2;
      break;
    case 1:
      minVolumes=kMinVoxelVolumesLevel3;
      break;
    default:
      minVolumes=10000;   // catch refinedDepth=3 and errors
      break;
  }

  if (refinedDepth<2)
  {
    std::size_t targetNo, noContainedDaughters, minNo, maxNo, replaceNo, i;
    G4double sliceWidth = (fmaxExtent-fminExtent)/maxNode;
    G4VoxelLimits newLimits;
    G4SmartVoxelNode* targetNode;
    G4SmartVoxelProxy* targetNodeProxy;
    G4SmartVoxelHeader* replaceHeader;
    G4SmartVoxelProxy* replaceHeaderProxy;
    G4VolumeNosVector* targetList;
    G4SmartVoxelProxy* lastProxy;
      
    for (targetNo=0; targetNo<maxNode; ++targetNo)
    {
      // Assume all slices are nodes (see preconditions)
      //
      targetNodeProxy = fslices[targetNo];
      targetNode = targetNodeProxy->GetNode();

      if (targetNode->GetNoContained() >= minVolumes)
      {
        noContainedDaughters = targetNode->GetNoContained();
        targetList = new G4VolumeNosVector();
        if (targetList == nullptr)
        {
          G4Exception("G4SmartVoxelHeader::RefineNodes()",
                      "GeomMgt0003", FatalException,
                      "Target volume node list allocation error.");
          return;
        }
        targetList->reserve(noContainedDaughters);
        for (i=0; i<noContainedDaughters; ++i)
        {
          targetList->push_back(targetNode->GetVolume((G4int)i));
        }
        minNo = targetNode->GetMinEquivalentSliceNo();
        maxNo = targetNode->GetMaxEquivalentSliceNo();

#ifdef G4GEOMETRY_VOXELDEBUG
        G4cout << "**** G4SmartVoxelHeader::RefineNodes" << G4endl
               << "     Refining nodes " << minNo 
               << " - " << maxNo << " inclusive" << G4endl;
#endif
        if (minNo > maxNo)    // Delete node and list to be replaced
        {                     // and avoid further action ...
          delete targetNode;
          delete targetList;
          return;
        }

        // Delete node proxies at start of collected sets of nodes/headers
        //
        lastProxy=nullptr;
        for (replaceNo=minNo; replaceNo<=maxNo; ++replaceNo)
        {
          if (lastProxy != fslices[replaceNo])
          {
            lastProxy=fslices[replaceNo];
            delete lastProxy;
          }
        }
        // Delete node to be replaced
        //
        delete targetNode;

        // Create new headers + proxies and replace in fslices
        //
        newLimits = pLimits;
        newLimits.AddLimit(faxis,fminExtent+sliceWidth*minNo,
                           fminExtent+sliceWidth*(maxNo+1));
        replaceHeader = new G4SmartVoxelHeader(pVolume,newLimits,
                                               targetList,(G4int)replaceNo);
        if (replaceHeader == nullptr)
        {
          G4Exception("G4SmartVoxelHeader::RefineNodes()", "GeomMgt0003",
                      FatalException, "Refined VoxelHeader allocation error.");
          return;
        }
        replaceHeader->SetMinEquivalentSliceNo((G4int)minNo);
        replaceHeader->SetMaxEquivalentSliceNo((G4int)maxNo);
        replaceHeaderProxy = new G4SmartVoxelProxy(replaceHeader);
        if (replaceHeaderProxy == nullptr)
        {
          G4Exception("G4SmartVoxelHeader::RefineNodes()", "GeomMgt0003",
                      FatalException, "Refined VoxelProxy allocation error.");
          return;
        }
        for (replaceNo=minNo; replaceNo<=maxNo; ++replaceNo)
        {
          fslices[replaceNo] = replaceHeaderProxy;
        }
        // Finished replacing current `equivalent' group
        //
        delete targetList;
        targetNo=maxNo;
      }
    }
  }
}

// ***************************************************************************
// Returns true if all slices have equal contents.
// Preconditions: all equal slices have been collected.
// Procedure:
// o checks all slice proxy pointers are equal
// o returns true if only one slice or all slice proxies pointers equal.
// ***************************************************************************
//
G4bool G4SmartVoxelHeader::AllSlicesEqual() const
{
  std::size_t noSlices = fslices.size();
  G4SmartVoxelProxy* refProxy;

  if (noSlices>1)
  {
    refProxy=fslices[0];
    for (std::size_t i=1; i<noSlices; ++i)
    {
      if (refProxy!=fslices[i])
      {
        return false;
      }
    }
  }
  return true;
}

// ***************************************************************************
// Streaming operator for debugging.
// ***************************************************************************
//
std::ostream& operator << (std::ostream& os, const G4SmartVoxelHeader& h)
{
  os << "Axis = " << G4int(h.faxis) << G4endl;
  G4SmartVoxelProxy *collectNode=nullptr, *collectHead=nullptr;
  std::size_t collectNodeNo = 0;
  std::size_t collectHeadNo = 0;
  std::size_t i, j;
  G4bool haveHeaders = false;

  for (i=0; i<h.fslices.size(); ++i)
  {
    os << "Slice #" << i << " = ";
    if (h.fslices[i]->IsNode())
    {
      if (h.fslices[i]!=collectNode)
      {
        os << "{";
        for (std::size_t k=0; k<h.fslices[i]->GetNode()->GetNoContained(); ++k)
        {
          os << " " << h.fslices[i]->GetNode()->GetVolume((G4int)k);
        }
        os << " }" << G4endl;
        collectNode = h.fslices[i];
        collectNodeNo = i;
      }
      else
      {
        os << "As slice #" << collectNodeNo << G4endl;
      }
    }
    else
    {
      haveHeaders=true;
      if (h.fslices[i] != collectHead)
      {
        os << "Header" << G4endl;
        collectHead = h.fslices[i];
        collectHeadNo = i;
      }
      else
      {
        os << "As slice #" << collectHeadNo << G4endl;
      }
    }
  }

  if (haveHeaders)
  {
    collectHead=nullptr;
    for (j=0; j<h.fslices.size(); ++j)
    {
      if (h.fslices[j]->IsHeader())
      {
        os << "Header at Slice #" << j << " = ";
        if (h.fslices[j] != collectHead)
        {
          os << G4endl 
             << (*(h.fslices[j]->GetHeader()));
          collectHead = h.fslices[j];
          collectHeadNo = j;
        }
        else
        {
          os << "As slice #" << collectHeadNo << G4endl;
        }
      }
    }
  }
  return os;
}
