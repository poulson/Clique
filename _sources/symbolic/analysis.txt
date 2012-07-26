Analysis
========
Once the elimination tree has been decided, it is necessary to determine 
the fill-in that will occur during factorization; that is to say, we must 
determine which degrees of freedom each (super)node of the elimination tree is 
coupled to after factorization. The set of vertices that a (super)node is 
coupled to is commonly referred to as its *structure*, and 
*symbolic factorization* consists of converting the original structure of each 
node in the elimination tree of a sparse matrix into the resulting structure
after factorization.

The following routine is called ``SymmetricAnalysis`` instead of 
``SymmetricSymbolicFactorization`` because, in addition to performing the 
symbolic factorization, it also computes and stores other data which is 
useful for numerical factorization and solves. 

.. note:: 
   Most users will not need to directly call the following routine, as it is 
   part of ``NestedDissection``. It will only need to be manually called in 
   cases where someone has manually computed the elimination tree of their 
   sparse matrix.

.. cpp:function:: void SymmetricAnalysis( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvIndices=true )
    
