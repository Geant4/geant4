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

/** An iterator that can be used to modify the tree. Most of the methods
 * provided by this class perform actions relative to the current position
 * of the iterator.
 */
class iterator : public tree_iterator<N, N&, N*>
{
  friend class tree;

  private:

    /* this typedef is very convenient and allows Java-like code */
    typedef tree_iterator<N, N&, N*> super;

  protected:

    iterator(tree_node<N>* current) : super(current) { }

  public:

    /** Creates an iterator that doesn't point to anything yet. */
    iterator(void) { }

    /** Copy constructor */
    iterator(const iterator& iter) : super(iter) { }

    /** Converts an iterator to a const_iterator.  This conversion
     * is one-way; const_iterator does not have a corresponding
     * conversion function to iterator.
     */
    operator const_iterator() const
    { return const_iterator(this->current); }

    typename super::reference
    operator * (void) const { return this->current->data; }

    typename super::pointer
    operator -> (void) const { return &(this->current->data); }

    /** Identical to the splice_after method except deletes the tree
     * after performing the splice.
     *
     * @return An iterator pointing to the root of the absorbed tree.
     */
    iterator absorb_after(tree<N>* t)
    {
      assert(t != NULL);
      assert(!t->empty());

      iterator i = splice_after(*t);
      delete t;
      return i;
    }

    /** Identical to the splice_back method except deletes the tree
     * after performing the splice.
     *
     * @return An iterator pointing to the root of the absorbed tree.
     */
    iterator absorb_back(tree<N>* t)
    {
#if !defined(TREE_NO_ERROR_CHECKING) && defined(TREE_EXCEPTIONS)
      if (t == NULL || t->empty())
      {
        throw typename tree<N>::null_tree_exception("in absorb_back");
      }
#else
      assert(t != NULL);
      assert(!t->empty());
#endif
 
      iterator i = splice_back(*t);
      delete t;
      return i;
    }

    /** Identical to the splice_before method except deletes the tree
     * after performing the splice.
     *
     * @return An iterator pointing to the root of the absorbed tree.
     */
    iterator absorb_before(tree<N>* t)
    {
      assert(t != NULL);
      assert(!t->empty());

      iterator i = splice_before(*t);
      delete t;
      return i;
    }

    /** Identical to the splice_front method except deletes the tree
     * after performing the splice.
     *
     * @return An iterator pointing to the root of the absorbed tree.
     */
    iterator absorb_front(tree<N>* t)
    {
      assert(t != NULL);
      assert(!t->empty());

      iterator i = splice_front(*t);
      delete t;
      return i;
    }

    sibling_iterator beginChildren(void) const
    { return sibling_iterator(iterator(this->current->first_child)); }

    /** Removes all children and grandchildren of this node
     */
    void clear(void);

    sibling_iterator endChildren(void) const
    { return sibling_iterator(iterator(NULL)); } 

    /** Removes this node and returns an iterator pointing
     * to the next node. */
    iterator erase(void);

    /** Inserts a new node before this node and returns an
     * iterator pointing to the new node. */
    iterator insert(const N& data);

    iterator insert(const tree<N>& t);

    /** Inserts a new node after this node and returns an
     * iterator pointing to the new node. */
    iterator insert_after(const N& data);

    iterator insert_after(const tree<N>& t); 

    iterator parent(void) const
    { return this->current->parent; }

    /** Removes the last child of this node.
     */
    void pop_back(void);

    /** Removes the first child of this node.
     */
    void pop_front(void);

    /** Adds a node to the end of this position's child list.
     *
     * @param data Data for the new node.
     *
     * @return An iterator pointing to the new node. 
     */
    iterator push_back (const N& data);

    /** Adds a copy of a tree to the end of this position's child list.
     *
     * @param t The tree to copy.
     *
     * @return An iterator pointing to the root of the inserted subtree.
     */
    iterator push_back (const tree<N>& t);

    /** Adds a node to the beginning of this position's child list.
     *
     * @param data Data for the new node.
     *
     * @return An iterator pointing to the new node. 
     */
    iterator push_front(const N& data);

    /** Adds a copy of a tree to the beginning of this position's child list.
     *
     * @param t The tree to copy.
     *
     * @return An iterator pointing to the root of the inserted subtree.
     */
    iterator push_front(const tree<N>& t);

    /** Replaces the data at this position.
     *
     * @param data The new data.
     *
     * @return The old data.
     */
    N replace(const N& data);

    /** Determines the size of the subtree below and including this position.
     *
     * @return The size of this subtree.
     */
    size_t size(void) const;

    /** Moves nodes from a tree (without copying) such that
     * the tree becomes the next sibling of this node.
     *
     * @param t The tree to move. It will be empty after this method returns.
     *
     * @return An iterator pointing to the root of the spliced tree.
     */
    iterator splice_after(tree<N>& t);

    /** Moves nodes from a tree (without copying) such that the tree becomes the last child of this node.
     *
     * @param t The tree to move. It will be empty after this method returns.
     *
     * @return An iterator pointing to the root of the spliced tree.
     */ 
    iterator splice_back(tree<N>& t);

    /** Moves nodes from a tree (without copying) such that the tree becomes the previous child of this node.
     *
     * @param t The tree to move. It will be empty after this method returns.
     *
     * @return An iterator pointing to the root of the spliced tree.
     */
    iterator splice_before(tree<N>& t);

    /** Moves nodes from a tree (without copying) such that the tree becomes the first child of this node.
     *
     * @param t The tree to move. It will be empty after this method returns.
     *
     * @return An iterator pointing to the root of the spliced tree.
     */
    iterator splice_front(tree<N>& t);
};
