//////////////////////////////////
/////  G4GAGTreeList.hh  ///////
/////   Satoshi Tanaka   ///////
//////////////////////////////////

#if !defined G4GAGTREELIST_HH
#define G4GAGTREELIST_HH

#include "g4std/iostream"	
#include "globals.hh"

/////////////////////////////////
/////  class G4GAGTreeNode  /////
/////////////////////////////////

template <class Type> class G4GAGTreeNode {
 public:
	Type item ;		// item
	G4GAGTreeNode<Type>* down ;	// address of neighboring node (downward)
	G4GAGTreeNode<Type>* up ;	// address of neighboring node (upward)
} ; // G4GAGTreeNode


///////////////////////////////////////
/////  class G4GAGTreeList<Type>  /////
///////////////////////////////////////

template <class Type> class G4GAGTreeList {
 private:
		//----- constants
	enum { ERROR = 0, NORMAL = 1 };

	G4GAGTreeNode<Type>*	head ;	// head of list
	G4GAGTreeNode<Type>*	tail ;	// tail of list
	G4GAGTreeNode<Type>*	cur ;	// current node

		//----- number of items in list
	int		nItem ; // number of stored items	

 public:
		//----- Constructor
		//.....  Initialize head and current node of list to NULL.
		//.....  Set name of list ( default name is "no_name")
	G4GAGTreeList( const char* name_ = "no_name" ) ;

		//----- Destructor
	~G4GAGTreeList() { this->Clear() ; }
	
		//----- Add an item to list and make it a new head.
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROOR = 0  (error)
	int	AddItem( const Type& item );	

		//----- Delete an item of current node
		//.....  Return: NORMAL = 1  if deleted.
		//.....          ERROR  = 0  if not deleted ( cur = NULL or nItem = 1).
		//.....  In case of NORMAL, cur moves to cur->down.
	int	DeleteItem() ;


		//----- Set current node to the head of list.
	void	ToHead()    { cur = head ; }

		//----- Set current node to the tail of list.
	void	ToTail()    { cur = tail ; }

		//----- Increment current node to downward direction by 1 step
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROOR = 0  (if cur = NULL)
	int	Downward();		

		//----- Increment current node to upward direction by 1 step
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROOR = 0  (if cur = NULL)
	int	Upward();		


		//----- Get current item
		//.....  The obtained item is returned to the argument.
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROR = 0   (error)
	int	GetItem( Type& item );

		//----- Get position of current node counting from tail
		//.....  Return: 0, 1, 2, ..... , nItem -1
		//.....  If cur == tail, 0 is returned.
		//.....  If cur == NULL, -1 is returned.
		//.....  Current node is not modified after calling this function
	int WhereIsCurrentNode() ;

		//----- Check if current node is NULL
		//.....  Return 1: if cur == NULL
		//.....         0: if cur != NULL
	int IsCurrentNodeNull() { return ( cur == NULL ) ; }

		//----- Get number of stored items
	int GetNumItem() const { return nItem ; }

		//----- Get number of stored items minus one
                //.....  Note: If the list is empty, -1 is returned
	int GetHeadIndex() const { return (nItem - 1) ; }

		//----- delete all items
		//..... All nodes are deleted.
		//..... head, tail, cur are reset to NULL.
	void Clear();

		//----- For use of STACK: Push
		//.....  Added new node becomes a new current node.
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROOR = 0  (error)
		//.....  Added new node becomes a new current node.
	int	Push( const Type& item ) { return AddItem( item ) ; }

		//----- For use of STACK
		//.....  Get head->item and then delete it.
		//.....  The obtained item is returned to the argument.
		//.....  Return: NORMAL = 1  (normal)
		//.....          ERROR = 0   (error)
		//.....  In case of NORMAL, cur moves to cur->down.
	int	Pop( Type& item ) ;

		//----- For use of STACK: pop
	void    Pop( void ) { this->ToHead(); this->DeleteItem(); }

} ; // class G4GAGTreeList



////////////////////////////////////////////////////////////////
/////  Inline Member functions of class "G4GAGTreeList"  ///////
////////////////////////////////////////////////////////////////

//----- 
template <class Type>
inline G4GAGTreeList<Type>::G4GAGTreeList( const char* name_ )
{
	cur = head = tail = NULL ;
	nItem = 0 ;

} // G4GAGTreeList<Type>::G4GAGTreeList()


//----- 
template <class Type>
inline int G4GAGTreeList<Type>::AddItem( const Type& item )
{
		//----- local
	int		status ;
	G4GAGTreeNode<Type>*	old_head = head ; // store old head

		//----- allocate dynamical memory for added node
	head = new G4GAGTreeNode<Type> ;		
	if( !head ){
			//----- failed to allocate dynamical memory
		status = ERROR ;

	} else {
			//----- increment number of stored items
		nItem++ ;				
	
			//----- set tail of list in the first push 
		if( nItem == 1 ) { tail = head ; }
	
			//----- set current node to added node
		cur = head ;
	
			//----- add an item
		head->item = item ;	// add an item
	
			//----- make links
			//.....          2          1
			//.....  NULL <----- head -----> old_head
			//.....                   <-----
			//.....                     3
		head->down = old_head ; // make link 1
		head->up   = NULL ;	// make link 2
		if( old_head != NULL ) {
			old_head->up = head ; // make link 3
		}
			//----- set status	
		status = NORMAL ;
	} // if--else

		//----- return status
	return status ;		

} 


//----- 
template <class Type>
inline  int G4GAGTreeList<Type>::DeleteItem()
{
		//----- local variables
	int status ;
	G4GAGTreeNode<Type>* cur_old ;	// backup of current node

		//----- ERROR if current node does not have an item
		//...... or no item exist in list
	if ( cur == NULL || nItem == 0 ) {
		status = ERROR ;
	} else {
			//----- Below cur != NULL, nItem >= 1, 
			//..... head != NULL, and tail != NULL
	
			//----- make backup of current node
		cur_old  = cur ;
	
			//----- remake links
		if( cur == head ) {
			cur  = cur->down ;
			head = cur ;
			if ( head != NULL ) {
				head->up = NULL ;
				delete cur_old ;  nItem-- ;
			} else {
				//----- In case that head = NULL:
				//.....  This case occurs in case that
				//.....  head = tail i.e. nItem = 1.
				//.....  At the end of this function,
				//.....  head = tail = cur = NULL, nItem = 0
				tail = NULL ;
				delete cur_old ; 
				nItem = 0 ;
				// nItem-- ;
			}
	
		} else if( cur == tail ){
			//----- Here, cur != head && cur == tail.
			//.....  So, tail != head.
			//.....  It means, at least 2 items exist in the list.
			//.....  So, tail->up is not NULL.
			tail = tail->up ;
			tail->down = NULL;
			cur = NULL ;
			delete cur_old ;  nItem-- ;
	
		} else {	 //----- In case that cur is neither head or tail.
				 //.....  head --- .... --- cur ---... --- tail.
			(cur->up)->down   = cur->down ;
			(cur->down)->up   = cur->up   ;
			cur = cur->down ;
			delete cur_old ;  nItem-- ;
		} // if-else_if-else
		
		status = NORMAL;

	} // if--else

		//----- return (successfully deleted)
	return status ;

} 


//----- 
template <class Type>
inline  int G4GAGTreeList<Type>::Downward()
{
	int status ;
	if( cur == NULL ) {
		status = ERROR ; // do nothing

	} else {
		cur = cur->down ; // move dawnward
		status = NORMAL ;
	}
	return status ;
} 


//----- 
template <class Type>
inline  int G4GAGTreeList<Type>::Upward()
{
	int status ;
	if( cur == NULL ) {
		status = ERROR ; // do nothing
	} else {
		cur = cur->up ; // move upward
		status = NORMAL ;
	}
	return status ;
} 


//----- 
template <class Type>
inline  int G4GAGTreeList<Type>::GetItem( Type& item )
{
		//----- return value
	int state ;

		//----- get an item of current node
	if ( cur == NULL ) {
		state = ERROR ;
	} else {
		state = NORMAL ;
		item = cur->item ;
	}

		//----- return
	return state ;
} 


//----- 
template <class Type>
inline  int G4GAGTreeList<Type>::WhereIsCurrentNode()
{
	G4GAGTreeNode<Type>* tmp = cur ;// copy current node
	int n = -1 ;			// counter

		//----- count how many steps there are until tail.
	while( tmp != NULL ) {	n++ ; tmp = tmp->down ;}

		//----- return steps
	return n ;

} 


//----- 
template <class Type>
inline  void G4GAGTreeList<Type>::Clear()
{
		//----- delete items
	cur = head ;
	while( cur != NULL )
	{
		G4GAGTreeNode<Type>* tmp = cur ;
		cur = cur->down ;
		delete tmp ;	// destructor called
	}

		//----- reset head and current node
	cur = head  = tail = NULL ;

		//----- reset number of stored items
	nItem = 0 ;

} 


//----- Pop()
template <class Type>
inline  int G4GAGTreeList<Type>::Pop( Type& item )
{
		//----- local variables
	int	status ;

		//----- get head->item
	ToHead() ;
	status = GetItem( item ) ;

		//----- delete head->item if it exists
		//.....  and make cur move to cur->down
	if( status ) { DeleteItem(); }

		//----- RETURN
	return status ;

} 



//////////////
#endif
////////////// end of list.h
