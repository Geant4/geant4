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
/*
 * G4Octree.cc
 *
 *  Created on: Feb 15, 2019
 *      Author: HoangTran
 */

#ifndef G4Octree_hh
#define G4Octree_hh 1
#include <array>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <iterator>
#include <iostream>
#include <typeinfo>
#include <list>
#include "G4ThreeVector.hh"
#include "G4DNABoundingBox.hh"
#include "G4Allocator.hh"

//using std::vector;
//using std::array;
//using namespace std;

const size_t max_per_node = 2;
const size_t max_depth = 100;

template <typename Iterator, class Extractor,typename Point = G4ThreeVector>
class G4Octree {
public:
    G4Octree();
    G4Octree(Iterator,Iterator);
    G4Octree(Iterator,Iterator, Extractor);

    using tree_type = G4Octree<Iterator,Extractor,Point>;

    //G4Octree(const tree_type& rhs);
    G4Octree(tree_type && rhs);
    void swap(tree_type& rhs);

    tree_type& operator = (tree_type rhs);
    tree_type& operator = (tree_type && rhs);

    ~G4Octree();

    size_t size() const;

    template <typename OutPutContainer>
    void radiusNeighbors(const Point& query, const G4double& radius,  OutPutContainer& resultIndices) const;

    void *operator new(size_t);
    void operator delete(void *);

private:
    enum NodeTypes
    {
        DEFAULT,
        LEAF,
        MAX_DEPTH_LEAF,
        INTERNAL
    };
    
    class Node;

    using NodeVector = std::vector<std::pair<Iterator,Point>>;
    using childNodeArray = std::array<Node*,8>;
    struct LeafValues
    {
        std::array<std::pair<Iterator,Point>,
        max_per_node> values_;
        size_t size_;
    };

    class Node
    {
    public:
        Node(const NodeVector& input_values);
        Node(const NodeVector& input_values,
             const G4DNABoundingBox& box,
             size_t current_depth);
        Node() = default;
        Node(const Node&) = delete;
        ~Node();
        template <typename OutPutContainer>
        G4bool radiusNeighbors(const Point& query, G4double radius,
                    OutPutContainer& resultIndices) const;
    private:
        void* fpValue;
        G4DNABoundingBox fBigVolume;
        NodeTypes fNodeType;

        void init_max_depth_leaf(const NodeVector& input_values);
        void init_leaf(const NodeVector& input_values);
        void init_internal(
                const NodeVector& input_values,
                size_t current_depth);
        struct InnerIterator
        {
            using wrapped_type = typename NodeVector::const_iterator;
            wrapped_type it__;
            InnerIterator(wrapped_type it):it__(it)
            {}
            Point operator*() const
            {
                return ((*it__).second);
            }
            InnerIterator& operator++()
            {
                ++it__;
                return *this;
            }
            InnerIterator operator++(G4int)
            {
                InnerIterator other = *this;
                ++it__;
                return other;
            }

            G4bool operator==(const InnerIterator& rhs) const
            {
                return this->it__ == rhs.it__;
            }

            G4bool operator!=(const InnerIterator& rhs) const
            {
                return !operator==(rhs);
            }
        };
    };
    private:
    Extractor functor_;
    Node* head_;
    size_t size_;
    G4ThreadLocalStatic G4Allocator<tree_type>* fgAllocator;
    
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include "G4Octree.icc"
#endif
