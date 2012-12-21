/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

template<typename T>
class SparseMatrix
{
public:
    // Construction and destruction
    SparseMatrix();
    SparseMatrix( int height );
    SparseMatrix( int height, int width );
    SparseMatrix( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    SparseMatrix( const DistSparseMatrix<T>& A );
    ~SparseMatrix();

    // High-level information
    int Height() const;
    int Width() const;
    const cliq::Graph& Graph() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numEntries );
    void Update( int row, int col, T value );
    int Capacity() const;

    // Data
    int Row( int entry ) const;
    int Col( int entry ) const;
    T Value( int entry ) const;
    int NumEntries() const;
    int EntryOffset( int row ) const;
    int NumConnections( int row ) const;

    // For modifying the size of the matrix
    void Empty();
    void ResizeTo( int height, int width );

    // For copying one matrix to another
    const SparseMatrix<T>& operator=( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    const SparseMatrix<T>& operator=( const DistSparseMatrix<T>& A );

    void Print( std::string msg ) const;

private:
    cliq::Graph graph_;
    std::vector<T> values_;

    template<typename U>
    struct Entry
    {
        int i, j;
        U value;
    };

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename U> friend class DistSparseMatrix;
};

} // namespace cliq
