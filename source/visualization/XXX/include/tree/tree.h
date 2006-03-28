/*
Copyright (c) 2003-2006, Troy Aaron Johnson
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  * Neither the name of Troy Aaron Johnson nor the names of any
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef TREE_H
#define TREE_H

#include <cassert>
#include <list>
#include <stdexcept>
#include <string>

/* master invariant checks are sprinkled throughout to check
   for catastrophic errors */
#if defined(TREE_NO_ERROR_CHECKING)
#define TREE_MASTER_INVARIANT
#else
#if defined(TREE_EXCEPTIONS)
#define TREE_MASTER_INVARIANT check_master_invariant()
#else
#define TREE_MASTER_INVARIANT assert(master != NULL); \
  assert(master->next_sibling != NULL)
#endif
#endif

/* experimental */
// #define TREE_POINTER_SPECIALIZATION

namespace taj  /* my initials */
{

/** Provides a generic n-ary acyclic tree container that behaves similarly
 * to and is mostly compatible with other standard C++ templates.  Various
 * iterators are provided, with pruning options.  Iterators are also used
 * to represent subtrees. 
 *
 * In the spirit of the standard template library, the BSD license makes
 * my tree class available for all to use.  I would appreciate receiving
 * patches to fix any bugs you may discover or suggestions on how to
 * improve the template.
 *
 * The STL list template was my initial inspiration for how to implement
 * the tree class, but I determined there was little benefit from
 * starting with that code.  This template class does not use any code
 * from the STL classes, but is designed to be compatible and have similar
 * method names for consistency.  I intend for this class to be very
 * efficient in terms of size and speed, just like the STL classes.
 * 
 * The std::list template uses many supporting classes that are found in
 * the std namespace or public in the std::list namespace.  This tree
 * template does not pollute the global namespace as much by keeping all
 * the supporting classes inside the main class.  If the user includes
 * the taj namespace, the only name that gets sucked into the enclosing
 * namespace is tree.  I feel this is a better design.  The following were
 * some motivating factors:
 *
 * 1) The tree_node class should be invisible to the user. The user wants to
 *    think of the N in tree<N> as the node type of the tree. In the standard
 *    list template, struct _List_node is similarly irrelevant to the user
 *    but unnecessarily appears in the std namespace.
 *
 * 2) The base template for the iterators should be a true class instead of a
 *    wide-open struct.  Furthermore, it should not be accessible to the user.
 *    A consequence of nesting the iterators inside the tree class is needing
 *    to make the tree their friend, but this is an example of how to use
 *    friends correctly.
 *
 * There is a lot of code here, so I have divided it into several header
 * files and code files.  The only header the user needs to include is tree.h.
 *
 *   tree.h                       - the main header file for the tree class
 *   tree_node.h                  - class that represents tree nodes
 *   tree_iterator.h              - base class for ALL iterators
 *   iterator.h                   - base class for ONLY non-const iterators
 *   const_iterator.h             - base class for ONLY const_iterators
 *   (const_)preorder_iterator.h  - iterators for preorder traversal
 *   (const_)postorder_iterator.h - iterators for postorder traversal
 *   (const_)bfs_iterator.h       - iterators for breadth-first traversal
 *   (const_)sibling_iterator.h   - iterators for a single tree level
 *
 * There is no "in-order" traversal because the tree is not necessarily
 * binary.  There are also .cc files which get included by this file, because
 * all the code for the template needs to be in tree.h.
 *
 * Note concerning GCC 3.4: Template code checks are much stricter in 3.4.
 * It is possible to have a template that compiles with GCC 3.3 but that does
 * not compile with 3.4.  The following GCC "bug" reports explain the
 * workarounds, which have been incorporated into this template class.
 *
 *   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=15552
 *   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=17649
 *
 * The same issue has found its way into the C++ FAQ
 *
 *  http://www.parashift.com/c++-faq-lite/templates.html#faq-35.12
 *  http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13
 *
 * A final note: C++ template code is weird and classes as complicated as
 * this one will give your compiler good exercise.  I have tested this
 * class with the Debian releases of GCC 3.3.5, GCC 3.4.4, GCC 4.0.3 and
 * I expect to use it with later versions.  I have no idea if it works
 * under other C++ compilers.  If there are changes that will allow it to
 * work under more compilers while still allowing it to work under GCC,
 * then you are encouraged to suggest them.
 *
 * @author Troy A. Johnson
 */
template <class N>
class tree
{
  private:

    #include "tree_node.h"
    #include "tree_iterator.h"

    /** The real root of the tree as opposed to the first node inserted by
     * the user. Also used as the end for preorder and postorder iterators
     * since they finish by "falling off" the root of the tree.  (The
     * breadth-first iterator uses NULL for its end because it falls off
     * the leaves.)  As far as the user is concerned, the first node is
     * master->next_sibling.  This corresponds to the _M_node in
     * _List_alloc_base for the standard list template.
     */
    tree_node<N>* master;

  public:

    /* avoid annoying forward reference problems */

    class iterator;
    class const_iterator;

    class sibling_iterator;
    class const_sibling_iterator;

    class preorder_iterator;
    class const_preorder_iterator;

    class postorder_iterator;
    class const_postorder_iterator;

    class bfs_iterator;
    class const_bfs_iterator;

    #include "iterator.h"
    #include "const_iterator.h"

    #include "sibling_iterator.h"
    #include "const_sibling_iterator.h"

    #include "preorder_iterator.h"
    #include "const_preorder_iterator.h"

    #include "postorder_iterator.h"
    #include "const_postorder_iterator.h"

    #include "bfs_iterator.h"
    #include "const_bfs_iterator.h"

    /* no dfs_iterator - use preorder_iterator instead */

    /** Thrown if tree exceptions are turned on and a NULL is encountered
     * in an incorrect place. */
    class null_tree_exception : public std::runtime_error
    {
      public:
        null_tree_exception(const std::string& s)
          : std::runtime_error("null_tree_exception " + s) { }
    };

  private:

#if !defined(TREE_NO_ERROR_CHECKING) && defined(TREE_EXCEPTIONS)
    void check_master_invariant(void)
    {
      if (master == NULL)
        throw new null_tree_exception("tree master is null");

      if (master->next_sibling == NULL)
        throw new null_tree_exception("tree master->next_sibling is null");
    }
#endif

    /** Initializes the master node.  Every tree constructor needs to
     * do this so it's a separate method.
     */
    void createMaster(void)
    {
      master = new tree_node<N>;
      master->parent = NULL;
      master->first_child  = master->last_child   = NULL;
      master->prev_sibling = master->next_sibling = master;
    }

    /** Replaces this tree with copies of all nodes below (and including)
     * the node pointed to by subtree.  Used by multiple constructors.
     *
     * @param subtree Location from which to begin copying.
     */
    void copy_subtree(const_iterator subtree);

    /** Determines the leftmost child of this tree.
     * Used as a starting point for postorder traversals.
     *
     * @return A pointer to the leftmost node or master if
     *    the tree is empty.
     */
    tree_node<N>* leftmostChild(void) const
    {
      TREE_MASTER_INVARIANT;

      tree_node<N>* p = master->next_sibling;

      while (p->first_child != NULL)
        p = p->first_child;

      return p;
    }

  public:

    /** Provides access to an iterator suitable for postorder traversal,
     * initially pointing to the leftmost child of this tree.
     *
     * @return An iterator for postorder traversal.
     */
    postorder_iterator       beginPost(void)
    { return postorder_iterator(leftmostChild()); }

    /** Provides access to an iterator suitable for postorder traversal,
     * initially pointing to the leftmost child of this tree.
     *
     * @return A constant iterator for postorder traversal.
     */
    const_postorder_iterator beginPost(void) const
    { return const_postorder_iterator(leftmostChild()); }
    
    /** Provides access to an iterator representing the end of a 
     * postorder traversal.
     *
     * @return The end of a postorder traversal.
     */
    postorder_iterator       endPost(void)
    { return postorder_iterator(master); }

    /** Provides access to an iterator representing the end of a 
     * constant postorder traversal.
     *
     * @return The end of a constant postorder traversal.
     */
    const_postorder_iterator endPost(void) const
    { return const_postorder_iterator(master); }

    /** Provides access to an iterator suitable for preorder traversal,
     * initially pointing to the root of this tree.
     *
     * @return An iterator for preorder traversal.
     */
    preorder_iterator        beginPre(void)
    { return preorder_iterator(master->next_sibling); }

    /** Provides access to an iterator suitable for preorder traversal,
     * initially pointing to the root of this tree.
     *
     * @return A constant iterator for preorder traversal.
     */
    const_preorder_iterator  beginPre(void) const
    { return const_preorder_iterator(master->next_sibling); }

    /** Provides access to an iterator representing the end of a 
     * preorder traversal.
     *
     * @return The end of a preorder traversal.
     */
    preorder_iterator        endPre(void)
    { return preorder_iterator(master); }

    /** Provides access to an iterator representing the end of a 
     * constant preorder traversal.
     *
     * @return The end of a constant preorder traversal.
     */
    const_preorder_iterator  endPre(void) const
    { return const_preorder_iterator(master); }

    /** Provides access to an iterator suitable for breadth-first traversal,
     * initially pointing to the root of this tree.
     *
     * @return An iterator for breadth-first traversal.
     */
    bfs_iterator             beginBfs(void)
    { return bfs_iterator(master->next_sibling); }

    /** Provides access to an iterator suitable for breadth-first traversal,
     * initially pointing to the root of this tree.
     *
     * @return A constant iterator for breadth-first traversal.
     */
    const_bfs_iterator       beginBfs(void) const
    { return const_bfs_iterator(master->next_sibling); }

    /** Provides access to an iterator representing the end of a
     * breadth-first traversal.
     *
     * @return The end of a breadth-first traversal.
     */
    bfs_iterator             endBfs(void)
    { return bfs_iterator(NULL); }

    /** Provides access to an iterator representing the end of a
     * constant breadth-first traversal.
     *
     * @return The end of a constant breadth-first traversal.
     */
    const_bfs_iterator       endBfs(void) const
    { return const_bfs_iterator(NULL); }
    
    /** Creates an empty tree.  A root node can be created with setRoot later.
     */
    tree(void)
    {
      createMaster();
    }

    /** Creates a single-node tree.
     *
     * @param root_data The data to use for the root of the tree.
     */
    explicit tree(const N& root_data)
    {
      createMaster();
      setRoot(root_data); 
    } 

    /** Copies an existing subtree.
     *
     * @param subtree The root of the original subtree
     *    from which to make the copy.
     */
    explicit tree(const_iterator subtree)
    {
      createMaster();
      copy_subtree(subtree);
    }
    
    /** Copies an existing tree.
     *
     * @param orig The original tree from which to make the copy.
     */
    tree(const tree<N>& orig)
    {
      createMaster();
      copy_subtree(orig.getRoot());
    }

    /** Creates a tree using copies of items of another container.
     * The range copied is [first, last).
     *
     * @param first First item to put in the tree.
     * @param last The item after the final item to put in the tree.
     * @param n Creates a n-ary tree.  Must be greater than zero.
     */
    template <class InputIterator>
    tree(InputIterator first, InputIterator last, unsigned int n);

    /** Copies an existing tree.
     *
     * @param orig The original tree from which to make the copy.
     *
     * @return A reference to this tree to be used in chained assignments.
     */ 
    const tree<N>& operator = (const tree<N>& orig)
    {
      copy_subtree(orig.getRoot());
      return *this;
    }

    /** Empties and destroys the tree.
     */
    virtual ~tree(void)
    {
      TREE_MASTER_INVARIANT;
      clear();
      delete master;
    }

    /** Empties the tree.  All nodes are deleted.  If the node type
     * is a pointer, it is the user's responsibility to delete
     * the data to which they point before clearing the tree.
     */
    void clear(void);

    /** Checks if the tree is empty.
     *
     * @return true if the tree is empty, false otherwise.
     */
    bool empty(void) const
    {
      TREE_MASTER_INVARIANT;
      return (master->next_sibling == master);
    }

    /** Provides access to the root of the tree.
     *
     * @return An iterator pointing at the root node.
     */
    iterator getRoot(void)
    {
      TREE_MASTER_INVARIANT;
      return iterator(master->next_sibling);
    }

    /** Provides access to the root of an unmodifiable tree.
     *
     * @return An iterator pointing at the root node.
     */
    const_iterator getRoot(void) const
    {
      TREE_MASTER_INVARIANT;
      return const_iterator(master->next_sibling);
    }

    /** Determines the height of the tree. This operation is
     * general and is linear in the size of the tree, not
     * the height.
     *
     * @return The height of the tree.
     */
    size_t height(void) const;

    /** Determines the height of the tree using
     * the leftmost grandchild of the root. This operation is
     * guaranteed to be linear in the height of the tree, but
     * may not be the true height for some trees.
     *
     * @return The height of the leftmost grandchild of the root.
     */
    size_t height_leftmost(void) const;

    /** Sets the data of the root node, or creates one if it was not set
     * when the tree was created.
     *
     * @return An iterator pointing to the root node.
     */
    iterator setRoot(const N& root_data)
    {
      if (empty())
      {
        tree_node<N>* root = new tree_node<N>(root_data);
        root->parent = NULL; 
        root->first_child = root->last_child = NULL;
        root->next_sibling = root->prev_sibling = master;

        master->prev_sibling = master->next_sibling = root;
      }
      else
        *getRoot() = root_data;

      return getRoot(); 
    }

    /** Determines the number of nodes in the tree.
     * Unfortunately this requires a full tree traversal.
     * This could be expensive, but it's far simpler than
     * making iterators aware of what tree they are modifying
     * and updating the current size of that tree. For simple
     * insertion and deletion that might not be difficult,
     * but for more complex operations it very well might be.
     * Doing it this way makes size a known expensive operation
     * instead of adding overhead to a large number of other calls.
     * If the user knows how many nodes they have added or deleted
     * since their last size call, then they can compute the current
     * size without calling size. It should be noted that the
     * standard library list template works the same way.
     *
     * @return The number of nodes in the tree.
     */
    size_t size(void) const;
};

#include "tree_iterator.cc"

#include "iterator.cc"

#include "const_iterator.cc"

/* keep this include last */
#include "tree.cc"

/* TODO - Pointer specialization

   This is experimental and may not work at all.
*/

#if defined(TREE_POINTER_SPECIALIZATION)
template class tree<void*>;

template <class N>
class tree<N*> : private tree<void*>
{
  public:

    /* do not prefix these forward declarations with tree<N*>:: */

    class iterator;
    class const_iterator;

    class sibling_iterator;
    class const_sibling_iterator;

    class preorder_iterator;
    class const_preorder_iterator;

    class postorder_iterator;
    class const_postorder_iterator;

    class bfs_iterator;
    class const_bfs_iterator;

    class iterator : private virtual tree<void*>::iterator
    {
      friend class tree<N*>;

      typedef class tree<N*>::iterator self;
      typedef class tree<void*>::iterator super;

      protected:

        iterator(const super& iter) : super(iter) { }

      public:

        iterator(void) { }

        iterator(const self& iter) : super(iter) { }

//        operator const_iterator()
//        { return static_cast<typename tree<void*>::const_iterator>(*this); }

        bool operator == (const self& iter) const
        { return super(*this) == super(iter); }

        bool operator != (const self& iter) const
        { return super(*this) != super(iter); }

        using super::clear;

        using super::depth;

        typename tree<N*>::sibling_iterator beginChildren(void) const
        { return self(super::beginChildren()); }

        typename tree<N*>::sibling_iterator endChildren(void) const
        { return self(super::endChildren()); }

        N*&
        operator * (void) const { return reinterpret_cast<N*&>(this->current->data); }

        N**
        operator -> (void) const { return reinterpret_cast<N**>(&(this->current->data)); }

        self absorb_back(tree<N*>* t)
        { return super::absorb_back(t); }

        self push_back(N* data)
        { return super::push_back(data); }

        N* replace(N* data)
        { return reinterpret_cast<N*>(super::replace(data)); }
    };

    class const_iterator : private virtual tree<void*>::const_iterator
    {
      friend class tree<N*>;

      typedef class tree<N*>::const_iterator self;
      typedef class tree<void*>::const_iterator super;

      protected:

        const_iterator(const super& iter) : super(iter) { }

      public:

        const_iterator(void) { }

        const_iterator(const self& iter) : super(iter) { }

        bool operator == (const self& iter) const
        { return super(*this) == super(iter); }

        bool operator != (const self& iter) const
        { return super(*this) != super(iter); }

        using super::depth;

        N* const &
        operator * (void) const { return reinterpret_cast<N* const &>(this->current->data); }

        N* const *
        operator -> (void) const { return reinterpret_cast<N* const *>(&(this->current->data)); }
    };

    class sibling_iterator : public tree<N*>::iterator, private tree<void*>::sibling_iterator
    {
      friend class tree<N*>;

      private:

        sibling_iterator(const tree<void*>::sibling_iterator& iter) : tree<N*>::iterator(iter) { }

      public:

        sibling_iterator(void) { }

        sibling_iterator(const tree<N*>::iterator& iter) : tree<N*>::iterator(iter) { }

        using tree<N*>::iterator::operator *;

        typename tree<N*>::sibling_iterator& operator ++ (void)
        { ++static_cast<tree<void*>::sibling_iterator>(*this); return *this; }
    };
#if 0
    class const_sibling_iterator : public const_iterator
    {
      friend class tree<N*>;

      private:

//        const_sibling_iterator(const tree<void*>::const_sibling_iterator& iter) : const_iterator(iter) { }

      public:

        const_sibling_iterator(const const_iterator& iter) : const_iterator(iter) { }
    };
#endif
    class preorder_iterator : public tree<N*>::iterator, private tree<void*>::preorder_iterator 
    {
      friend class tree<N*>;

      private:

        preorder_iterator(const tree<void*>::preorder_iterator& iter) : tree<N*>::iterator(iter) { }

      public:

        preorder_iterator(const tree<N*>::iterator& iter) : tree<N*>::iterator(iter) { }

        using tree<N*>::iterator::operator ==;
        using tree<N*>::iterator::operator !=;
        using tree<N*>::iterator::operator *;

        typename tree<N*>::preorder_iterator& operator ++ (void)
        { ++static_cast<tree<void*>::preorder_iterator>(*this); return *this; }

        using tree<N*>::iterator::replace;
    };

    class const_preorder_iterator : public tree<N*>::const_iterator, private tree<void*>::const_preorder_iterator
    {
      friend class tree<N*>;

      private:

        const_preorder_iterator(const tree<void*>::const_preorder_iterator& iter) : tree<void*>::const_preorder_iterator(iter) { }

      public:

        const_preorder_iterator(const tree<N*>::const_iterator& iter) : tree<N*>::const_iterator(iter) { }

        using tree<N*>::const_iterator::operator *;

        typename tree<N*>::const_preorder_iterator& operator ++ (void)
        { ++static_cast<tree<void*>::const_preorder_iterator>(*this); return *this; }
    };
#if 0
    class postorder_iterator : private tree<void*>::postorder_iterator
    {
      friend class tree<N*>;

      private:

//        postorder_iterator(const tree<void*>::postorder_iterator& iter) : iterator(iter) { }

      public:

        postorder_iterator(const iterator& iter) : tree<void*>::postorder_iterator(iter) { }
    };

    class const_postorder_iterator : private tree<void*>::const_postorder_iterator
    {
      friend class tree<N*>;

      private:

//        const_postorder_iterator(const tree<void*>::const_postorder_iterator& iter) : const_iterator(iter) { }

      public:

        const_postorder_iterator(const const_iterator& iter) : tree<void*>::const_postorder_iterator(iter) { }
    };

    class bfs_iterator : private tree<void*>::bfs_iterator
    {
      friend class tree<N*>;

      private:

//        bfs_iterator(const tree<void*>::bfs_iterator& iter) : iterator(iter) { }

      public:

        bfs_iterator(const iterator& iter) : tree<void*>::bfs_iterator(iter) { }
    };

    class const_bfs_iterator : private tree<void*>::const_bfs_iterator
    {
      friend class tree<N*>;

      private:

//        const_bfs_iterator(const tree<void*>::const_bfs_iterator& iter) : const_iterator(iter) { }

      public:

        const_bfs_iterator(const const_iterator& iter) : tree<void*>::const_bfs_iterator(iter) { }
    };
#endif
    tree(void) : tree<void*>::tree()
    {
    }

    /** Creates a single-node tree.
     *
     * @param root_data The data to use for the root of the tree.
     */
    explicit tree(N* root_data) : tree<void*>::tree(root_data)
    {
    }

    explicit tree(const_iterator subtree) : tree<void*>::tree(subtree)
    {
    }

    tree(const tree<N*>& orig) : tree<void*>::tree(orig)
    {
    }

    preorder_iterator        beginPre(void)
    { return tree<void*>::beginPre(); }

    const_preorder_iterator  beginPre(void) const
    { return tree<void*>::beginPre(); }

    preorder_iterator        endPre(void)
    { return tree<void*>::endPre(); }

    const_preorder_iterator  endPre(void) const
    { return tree<void*>::endPre(); }

    iterator getRoot(void)
    { return tree<void*>::getRoot(); }

    const_iterator getRoot(void) const
    { return tree<void*>::getRoot(); }

    iterator setRoot(N* root_data)
    { return tree<void*>::setRoot(root_data); }
};

#if 0
template class tree<void*>;
template class tree<void*>::tree_node<void*>;

template <class N>
class tree<N*> : public tree<void*>
{
  public:

//    class const_iterator;

  private:

    template <class T>
    class tree_node<T*> : public tree_node<void*>
    {
      public:

        tree_node(void) : tree_node<void*>() { }

        tree_node(T* data) : tree_node<void*>(data) { }
    }; 
/*
    template <class T, class R, class P>
    class tree_iterator<T*, R*, P*> : public tree_iterator<void*, void*&, void**>
    {

    };
*/
  public:

    tree(void) : tree<void*>()
    {
      std::cout << "specialized tree default constructor" << std::endl;
    }

    explicit tree(N* root_data) : tree<void*>(root_data)
    {
    }
};
/*
template<class N>
class tree<N*>::iterator : public tree<N>::template tree_iterator<N*, N&, N*>
{

};
*/

template <class N>
class tree<N*>::const_iterator : public tree<void*>::iterator
{
  private:

    /* this typedef is very convenient and allows Java-like code */
    typedef typename tree<N*>::template tree_iterator<N*, const N*&, const N**> super;

  public:

    typename super::reference
    operator * (void) const { return this->current->data; }
};

//template <class N*>
//typename tree<N*>::template tree_iterator<N*, const N*&, const N**>::reference
//tree<N*>::template const_iterator::operator * (void) const
//{ return this->current->data; }

#endif
#endif

} /* namespace taj */

#undef TREE_MASTER_INVARIANT

#endif
