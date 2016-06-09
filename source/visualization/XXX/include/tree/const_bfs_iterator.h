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

class const_bfs_iterator : public virtual const_iterator
{
  friend class tree;

  private:

    std::list<tree_node<N>*> queue;

    explicit const_bfs_iterator(tree_node<N>* current) : const_iterator(current) { }

  public:

    const_bfs_iterator(void) { }
    const_bfs_iterator(const const_iterator& iter) : const_iterator(iter) { }

    const_bfs_iterator endBfs(void) const { return const_bfs_iterator(NULL); }

    const_bfs_iterator& operator ++ (void)
    {
      /* this line necessary for g++ 3.4 onward */
      tree_node<N>*& current = this->current;

      assert(current != NULL);

      for (tree_node<N>* p = current->first_child;
           p != NULL; p = p->next_sibling)
      {
        queue.push_front(p);
      }

      if (queue.empty())
      {
        current = NULL;
      }
      else
      {
        current = queue.back();
        queue.pop_back();
      }

      return *this;
    }
};

