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

/** An iterator that can be used to inspect the tree but not modify it.
 */
class const_iterator : public tree_iterator<N, const N&, const N*>
{
  friend class tree;

  /* allows iterator to call the protected constructor below */
  friend tree<N>::iterator::operator const_iterator() const;

  private:

    /* this typedef is very convenient and allows Java-like code */
    typedef tree_iterator<N, const N&, const N*> super;

  protected:

    const_iterator(tree_node<N>* current) : super(current) { }

  public:

    const_iterator(void) { }
    const_iterator(const const_iterator& iter) : super(iter) { }

    typename super::reference
    operator * (void) const { return this->current->data; }

    typename super::pointer
    operator -> (void) const { return &(this->current->data); }

    const_sibling_iterator beginChildren(void) const
    { return const_sibling_iterator(const_iterator(this->current->first_child)); }

    const_sibling_iterator endChildren(void) const
    { return const_sibling_iterator(const_iterator(NULL)); }

    const_iterator parent(void) const
    { return this->current->parent; }

    /** Determines the size of the subtree below and including this position.
     *
     * @return The size of this subtree.
     */
    size_t size(void) const;
};
