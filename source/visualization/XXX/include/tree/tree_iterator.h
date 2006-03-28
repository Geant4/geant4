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

/** The base template for all iterators.
 */
template <class T, class R, class P>
class tree_iterator
{
  friend class tree;

  protected:

    /* The node pointed to by this iterator. */
    tree_node<N>* current; 

    /* This is the depth of an iterator relative to its starting position.
       The starting depth is 0, the depth of the parent is -1, the depth
       of any child is 1, the depth of any grandchild is 2, etc. */
    long iter_depth;

    /** Allows derived classes to create an iterator
     * from a pointer to a node.  This constructor
     * is not explicit because the tree_iterator class
     * contains the == and != operators, which are used
     * throughout the template to compare iterators
     * to NULL.  The tree_iterator class is not visible
     * to the user, so this implicit constructor is
     * entirely for my convenience when coding the template.
     *
     * @param current Position for the new iterator
     */ 
    tree_iterator(tree_node<N>* current) : current(current), iter_depth(0L) { }

    /** Destroys the node pointed to by this iterator.
     */
    void destroy(void)
    {
      delete current;
      current = NULL;
      iter_depth = 0L;
    }

    /* increment and decrement methods used by both the iterator
       and the const_iterator versions of various iterators */
  
    static inline void post_inc(tree_node<N>*& current, long& depth);
    static inline void post_dec(tree_node<N>*& current, long& depth);

    static inline void pre_inc(tree_node<N>*& current, long& depth);
    static inline void pre_dec(tree_node<N>*& current, long& depth);
 
    static inline void sib_inc(tree_node<N>*& current);
    static inline void sib_dec(tree_node<N>*& current);
 
  public:

    /* The following typedefs are required for compatibility
     * with standard containers. */
    
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef T                               value_type;
    typedef R                               reference;
    typedef P                               pointer;
    typedef size_t                          size_type;
    typedef ptrdiff_t                       difference_type;

    tree_iterator(void) : current(NULL), iter_depth(0L) { }

    tree_iterator(const tree_iterator& iter)
      : current(iter.current), iter_depth(iter.iter_depth) { }

    /* note that depth is not relevant when comparing iterators */

    bool operator == (const tree_iterator& iter) const
    { return this->current == iter.current; }

    bool operator != (const tree_iterator& iter) const
    { return this->current != iter.current; }

    long depth(void) const { return iter_depth; }

    bool leaf(void) const { return this->first_child == NULL; }
};
