//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SmartVoxelHeader.cc,v 1.18 2002-05-15 10:06:33 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4SmartVoxelHeader
//
// Implementation
//
// Define G4GEOMETRY_VOXELDEBUG for debugging information on G4cout
//
// History:
// 29.04.02 Use 3D voxelisation for non consuming replication - G.C.
// 18.04.01 Migrated to STL vector - G.C.
// 12.02.99 Introduction of new quality/smartless: max for (slices/candid) S.G.
// 11.02.99 Voxels at lower levels are now built for collapsed slices      S.G.
// 21.07.95 Full implementation, supporting non divided physical volumes
// 14.07.95 Initial version - stubb definitions only
// ***************************************************************************

#include "G4SmartVoxelHeader.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"

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
  G4int nDaughters = pVolume->GetNoDaughters();
  G4VoxelLimits limits;   // Create `unlimited' limits object

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
  for (size_t i=0;i<pCandidates->size();i++)
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
  G4int node, proxy, maxNode=fslices.size();
  G4SmartVoxelProxy *lastProxy=0;
  G4SmartVoxelNode *dyingNode, *lastNode=0;
  G4SmartVoxelHeader *dyingHeader, *lastHeader=0;

  for (node=0; node<maxNode; node++)
  {
    if (fslices[node]->IsHeader())
    {
      dyingHeader = fslices[node]->GetHeader();
      if (lastHeader!=dyingHeader)
      {
        lastHeader = dyingHeader;
        lastNode = 0;
        delete dyingHeader;
      }
    }
      else
    {
      dyingNode = fslices[node]->GetNode();
      if (dyingNode!=lastNode)
      {
        lastNode=dyingNode;
        lastHeader=0;
        delete dyingNode;
      }
    }
  }
  // Delete proxies
  //
  for (proxy=0; proxy<maxNode; proxy++)
  {
    if (fslices[proxy]!=lastProxy)
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
    G4int node, maxNode;
    G4SmartVoxelProxy *leftProxy, *rightProxy;
    G4SmartVoxelHeader *leftHeader, *rightHeader;
    G4SmartVoxelNode *leftNode, *rightNode;

    maxNode=GetNoSlices();
    for (node=0; node<maxNode; node++)
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
          if (!(*leftHeader==*rightHeader))
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
          if (!(*leftNode==*rightNode))
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
  G4int nDaughters = pVolume->GetNoDaughters();

  G4VolumeNosVector targetList;
  targetList.reserve(nDaughters);
  for (G4int i=0; i<nDaughters; i++)
  {
    targetList.push_back(i);
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
  G4VPhysicalVolume *pDaughter=0;

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
    pDaughter=pVolume->GetDaughter(0);
    pDaughter->GetReplicationData(axis,nReplicas,width,offset,consuming);
    fparamAxis = axis;
    if ( consuming==false )
    {
      G4VoxelLimits limits;   // Create `unlimited' limits object
      G4VolumeNosVector targetList;
      targetList.reserve(nReplicas);
      for (G4int i=0; i<nReplicas; i++)
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
          G4cout << "ERROR - Illegal axis !" << G4endl;
          G4Exception("ERROR - G4SmartVoxelHeader::BuildReplicaVoxels");
          break;
      }  
      faxis = axis;   // Set axis
      BuildConsumedNodes(nReplicas);
      if ( (axis==kXAxis) || (axis==kYAxis) || (axis==kZAxis) )
      {
        // Sanity check on extent
        //
        G4double min, max;
        G4VoxelLimits limits;
        G4AffineTransform origin;
        pVolume->GetSolid()->CalculateExtent(axis, limits, origin, min, max);
        if ( (fabs((min-fminExtent)/fminExtent) +
              fabs((max-fmaxExtent)/fmaxExtent)) > 0.05)
        {
          G4cout << "ERROR - Replicated geometry, logical volume: "
                 << pVolume->GetName() << G4endl;
          G4Exception("ERROR - G4SmartVoxelHeader::BuildReplicaVoxels");
        }
      }
    }
  }
  else
  {
    G4cout << "ERROR - There must be a single replicated volume !" << G4endl;
    G4Exception("ERROR - G4SmartVoxelHeader::BuildReplicaVoxels");
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
  G4SmartVoxelNode *pNode;
  G4SmartVoxelProxy *pProxyNode;

  // Create and fill nodes in temporary G4NodeVector (on stack)
  //
  G4NodeVector nodeList;
  nodeList.reserve(nReplicas);
  for (nNode=0; nNode<nReplicas; nNode++)
  {
    pNode=new G4SmartVoxelNode(nNode);
    if (!pNode)
    {
      G4cout << "ERROR - Node allocation failed." << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildConsumedNodes");
    }
    nodeList.push_back(pNode);
  }
  for (nVol=0; nVol<nReplicas; nVol++)
  {
    nodeList[nVol]->Insert(nVol);   // Insert replication of number
  }                                 // identical to voxel number

  // Create & fill proxy List `in place' by modifying instance data fslices
  //
  fslices.clear();
  for (nNode=0; nNode<nReplicas; nNode++)
  {
    pProxyNode = new G4SmartVoxelProxy(nodeList[nNode]);
    if (!pProxyNode)
    {
      G4cout << "ERROR - Proxy Node allocation failed." << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildConsumedNodes");
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

  G4ProxyVector *pGoodSlices=0, *pTestSlices, *tmpSlices;
  G4double goodSliceScore=kInfinity, testSliceScore;
  EAxis goodSliceAxis = kXAxis;
  EAxis testAxis      = kXAxis;
  G4int node, maxNode, iaxis;
  G4VoxelLimits noLimits;

  // Try all non-limited cartesian axes
  //
  for (iaxis=0; iaxis<3; iaxis++)
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
        for (node=0; node<maxNode; node++)
        {
          delete (*pTestSlices)[node]->GetNode();
        }
        G4SmartVoxelProxy* tmpProx;
        while (pTestSlices->size()>0)
        {
          tmpProx = pTestSlices->back();
          pTestSlices->pop_back();
          for (G4ProxyVector::iterator i=pTestSlices->begin();
                                       i!=pTestSlices->end(); i++)
          {
            if (*i==tmpProx)
            {
              pTestSlices->erase(i); i--;
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
    G4cout << "ERROR - Illegal limits: already 3 dimensions of limits !" << G4endl;
    G4Exception("ERROR - G4SmartVoxelHeader::BuildVoxelsWithinLimits");
  }

  // 
  // We have selected pGoodSlices, with a score testSliceScore
  //

  // Store chosen axis, slice ptr
  //
  fslices=*pGoodSlices; // Set slice information, copy ptrs in collection
  delete pGoodSlices;   // Destroy slices vector, but not contained
                        // proxies or nodes
  faxis=goodSliceAxis;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << G4endl << "     Selected axis = " << faxis << G4endl;
  for (size_t islice=0; islice<fslices.size(); islice++)
  {
    G4cout << "     Node #" << islice << " = {";
    for (G4int j=0; j<fslices[islice]->GetNode()->GetNoContained(); j++)
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
  CollectEquivalentNodes();     // Collect common nodes
  RefineNodes(pVolume,pLimits); // Refine nodes creating headers

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
  G4int sliceNo, minNo, maxNo, equivNo;
  G4int maxNode = fslices.size();
  G4SmartVoxelNode *startNode, *sampleNode;
  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
  {
    minNo = sliceNo;

    // Get first node (see preconditions - will throw exception if a header)
    //
    startNode = fslices[minNo]->GetNode();

    // Find max equivalent
    //
    for (equivNo=minNo+1; equivNo<maxNode; equivNo++)
    {
      sampleNode = fslices[equivNo]->GetNode();
      if (!((*startNode) == (*sampleNode))) { break; }
    }
    maxNo = equivNo-1;
    if (maxNo != minNo)
    {
      // Set min and max nos
      //
      for (equivNo=minNo; equivNo<=maxNo; equivNo++)
      {
        sampleNode = fslices[equivNo]->GetNode();
        sampleNode->SetMinEquivalentSliceNo(minNo);
        sampleNode->SetMaxEquivalentSliceNo(maxNo);
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
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode=fslices.size();
  G4SmartVoxelNode *equivNode;
  G4SmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
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
      for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
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
  G4int sliceNo, maxNo, equivNo;
  G4int maxNode = fslices.size();
  G4SmartVoxelHeader *equivHeader, *sampleHeader;
  G4SmartVoxelProxy *equivProxy;

  for (sliceNo=0; sliceNo<maxNode; sliceNo++)
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
        for (equivNo=sliceNo+1; equivNo<=maxNo; equivNo++)
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
  G4double motherMinExtent, motherMaxExtent, targetMinExtent, targetMaxExtent;
  G4VPhysicalVolume *pDaughter=0;
  G4VPVParameterisation *pParam=0;
  G4VSolid *targetSolid;
  G4AffineTransform targetTransform;
  G4bool replicated;
  G4int nCandidates = pCandidates->size();
  G4int nVol, nNode, targetVolNo;
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

  if ( (pVolume->GetNoDaughters()==1)
    && (pVolume->GetDaughter(0)->IsReplicated()==true) )
  {
    // Replication data not required: only parameterisation object 
    // and volume no. List used
    //
    pDaughter = pVolume->GetDaughter(0);
    pParam = pDaughter->GetParameterisation();
    if (!pParam)
    {
      G4cout << "PANIC! Replicated volume with no parameterisation object !"
             << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
    }

    // Setup volume, preserving current mother link
    //
    pDaughter->Setup(pDaughter->GetMother());
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
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    targetVolNo=(*pCandidates)[nVol];
    if (replicated == false)
    {
      pDaughter=pVolume->GetDaughter(targetVolNo);

      // Setup volume, preserving current mother link
      //
      pDaughter->Setup(pDaughter->GetMother());
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
      targetSolid = pParam->ComputeSolid(targetVolNo,pDaughter);

      // Setup solid
      //
      targetSolid->ComputeDimensions(pParam,targetVolNo,pDaughter);

      // Setup transform
      //
      pParam->ComputeTransformation(targetVolNo,pDaughter);
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

    // Check not entirely outside mother when processing toplevel nodes
    //
    if ( (!pLimits.IsLimited()) && ((targetMaxExtent<=motherMinExtent)
                                  ||(targetMinExtent>=motherMaxExtent)) )
    {
      G4cout << "PANIC! Daughter physical volume "
             << pDaughter->GetName() << G4endl
	     << "is entirely outside mother logical volume "
	     << pVolume->GetName() << " !!" << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
    }

#ifdef G4GEOMETRY_VOXELDEBUG
    // Check for straddling volumes when debugging.
    // If a volume is >kStraddlePercent percent over the mother
    // boundary, print a warning.
    //
    if (!pLimits.IsLimited())
    {
      G4double width;
      G4int kStraddlePercent=5;
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
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    currentWidth = maxExtents[nVol]-minExtents[nVol];
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
  G4int    noNodes = ((noNodesSmart-noNodesExactI)>=0.5)
                     ? noNodesExactI+1 : noNodesExactI;
  if( noNodes == 0 ) { noNodes=1; }

#ifdef G4GEOMETRY_VOXELDEBUG
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
  G4NodeVector* nodeList = new G4NodeVector();
  nodeList->reserve(noNodes);
  if (!nodeList)
  {
    G4cout << "ERROR - NodeList allocation failed." << G4endl;
    G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
  }
  for (nNode=0; nNode<noNodes; nNode++)
  {
    G4SmartVoxelNode *pNode;
    pNode = new G4SmartVoxelNode(nNode);
    if (!pNode)
    {
      G4cout << "ERROR - Node allocation failed." << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
    }
    nodeList->push_back(pNode);
  }

  // All nodes created (empty)

  // Fill nodes: Step through extent lists
  //
  for (nVol=0; nVol<nCandidates; nVol++)
  {
    G4int nodeNo, minContainingNode, maxContainingNode;
    minContainingNode = G4int((minExtents[nVol]-motherMinExtent)/nodeWidth);
    maxContainingNode = G4int((maxExtents[nVol]-motherMinExtent)/nodeWidth);

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
        minContainingNode=0;
      }
      for (nodeNo=minContainingNode; nodeNo<=maxContainingNode; nodeNo++)
      {
        (*nodeList)[nodeNo]->Insert((*pCandidates)[nVol]);
      }
    }
  }

  // All nodes filled

  // Create proxy List : caller has deletion responsibility
  // (but we must delete nodeList *itself* - not the contents)
  //
  G4ProxyVector* proxyList = new G4ProxyVector();
  proxyList->reserve(noNodes);
  if (!proxyList)
  {
    G4cout << "ERROR - Proxy List allocation failed." << G4endl;
    G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
  }
  //
  // Fill proxy List
  //
  for (nNode=0; nNode<noNodes; nNode++)
  {
    G4SmartVoxelProxy* pProxyNode = new G4SmartVoxelProxy((*nodeList)[nNode]);
    if (!pProxyNode)
    {
      G4cout << "ERROR - Proxy Node allocation failed." << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::BuildNodes");
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
  G4int nNodes = pSlice->size();
  G4int noContained, maxContained=0, sumContained=0, sumNonEmptyNodes=0;
  G4SmartVoxelNode *node;

  for (G4int i=0; i<nNodes; i++)
  {
    if ((*pSlice)[i]->IsNode())
    {
      // Definitely a node. Add info to running totals
      //
      node = (*pSlice)[i]->GetNode();
      noContained = node->GetNoContained();
      if (noContained)
      {
        sumNonEmptyNodes++;
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
      G4cout << "ERROR - Not defined for sliced volumes." << G4endl;
      G4Exception("ERROR - G4SmartVoxelHeader::CalculateQuality");
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
  G4int refinedDepth=0, minVolumes;
  G4int maxNode = fslices.size();

  if (pLimits.IsXLimited()) 
  {
    refinedDepth++;
  }
  if (pLimits.IsYLimited()) 
  {
    refinedDepth++;
  }
  if (pLimits.IsZLimited()) 
  {
    refinedDepth++;
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
    G4int targetNo, noContainedDaughters, minNo, maxNo, replaceNo, i;
    G4double sliceWidth = (fmaxExtent-fminExtent)/maxNode;
    G4VoxelLimits newLimits;
    G4SmartVoxelNode* targetNode;
    G4SmartVoxelProxy* targetNodeProxy;
    G4SmartVoxelHeader* replaceHeader;
    G4SmartVoxelProxy* replaceHeaderProxy;
    G4VolumeNosVector* targetList;
    G4SmartVoxelProxy* lastProxy;
      
    for (targetNo=0; targetNo<maxNode; targetNo++)
    {
      // Assume all slices are nodes (see preconditions)
      //
      targetNodeProxy = fslices[targetNo];
      targetNode = targetNodeProxy->GetNode();

      if (targetNode->GetNoContained() >= minVolumes)
      {
        noContainedDaughters = targetNode->GetNoContained();
        targetList = new G4VolumeNosVector();
        targetList->reserve(noContainedDaughters);
        if (!targetList)
        {
          G4cout << "ERROR - Target volume no List new failed." << G4endl;
          G4Exception("ERROR - G4SmartVoxelHeader::RefineNodes");
        }
        for (i=0; i<noContainedDaughters; i++)
        {
          targetList->push_back(targetNode->GetVolume(i));
        }
        minNo = targetNode->GetMinEquivalentSliceNo();
        maxNo = targetNode->GetMaxEquivalentSliceNo();

#ifdef G4GEOMETRY_VOXELDEBUG
        G4cout << "**** G4SmartVoxelHeader::RefineNodes" << G4endl
               << "     Refining nodes " << minNo 
               << " - " << maxNo << " inclusive" << G4endl;
#endif
        // Delete node proxies at start of collected sets of nodes/headers
        //
        lastProxy=0;
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
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
                                               targetList,replaceNo);
        if (!replaceHeader)
        {
          G4cout << "ERROR - Refined VoxelHeader new failed." << G4endl;
          G4Exception("ERROR - G4SmartVoxelHeader::RefineNodes");
        }
        replaceHeader->SetMinEquivalentSliceNo(minNo);
        replaceHeader->SetMaxEquivalentSliceNo(maxNo);
        replaceHeaderProxy = new G4SmartVoxelProxy(replaceHeader);
        if (!replaceHeader)
        {
          G4cout << "ERROR - Refined VoxelProxy new failed." << G4endl;
          G4Exception("ERROR - G4SmartVoxelHeader::RefineNodes");
        }
        for (replaceNo=minNo; replaceNo<=maxNo; replaceNo++)
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
  G4int noSlices = fslices.size();
  G4SmartVoxelProxy* refProxy;

  if (noSlices>1)
  {
    refProxy=fslices[0];
    for (G4int i=1; i<noSlices; i++)
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
G4std::ostream& operator << (G4std::ostream& s, const G4SmartVoxelHeader& h)
{
  s << "Axis = " << G4int(h.faxis) << G4endl;
  G4SmartVoxelProxy *collectNode=0, *collectHead=0;
  G4int collectNodeNo=0;
  G4int collectHeadNo=0;
  size_t i, j;
  G4bool haveHeaders=false;

  for (i=0; i<h.fslices.size(); i++)
  {
    s << "Slice #" << i << " = ";
    if (h.fslices[i]->IsNode())
    {
      if (h.fslices[i]!=collectNode)
      {
        s << "{";
        for (G4int j=0; j<h.fslices[i]->GetNode()->GetNoContained(); j++)
        {
          s << " " << h.fslices[i]->GetNode()->GetVolume(j);
	}
        s << " }" << G4endl;
        collectNode = h.fslices[i];
        collectNodeNo = i;
      }
      else
      {
        s << "As slice #" << collectNodeNo << G4endl;
      }
    }
    else
    {
      haveHeaders=true;
      if (h.fslices[i] != collectHead)
      {
        s << "Header" << G4endl;
        collectHead = h.fslices[i];
        collectHeadNo = i;
      }
      else
      {
        s << "As slice #" << collectHeadNo << G4endl;
      }
    }
  }

  if (haveHeaders)
  {
    collectHead=0;
    for (j=0; j<h.fslices.size(); j++)
    {
      if (h.fslices[j]->IsHeader())
      {
        s << "Header at Slice #" << j << " = ";
        if (h.fslices[j] != collectHead)
        {
          s << G4endl 
            << (*(h.fslices[j]->GetHeader()));
          collectHead = h.fslices[j];
          collectHeadNo = j;
        }
        else
        {
          s << "As slice #" << collectHeadNo << G4endl;
        }
      }
    }
  }
  return s;
}
