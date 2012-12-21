/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename T>
class DistVector
{
public:
    // Constructors and destructors
    DistVector();
    DistVector( mpi::Comm comm );
    DistVector( int height, mpi::Comm comm );
    DistVector( int height, T* buffer, mpi::Comm comm );
    DistVector( int height, const T* buffer, mpi::Comm comm );
    // TODO: Constructor for building from a DistVector
    ~DistVector();

    // High-level information
    int Height() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalRow() const;
    int LocalHeight() const;

    // Local data
    T GetLocal( int localRow ) const;
    void SetLocal( int localRow, T value );
    void UpdateLocal( int localRow, T value );
    Matrix<T>& LocalVector();
    const Matrix<T>& LocalVector() const;

    // For modifying the size of the vector
    void Empty();
    void ResizeTo( int height );

    // Assignment
    const DistVector<T>& operator=( const DistVector<T>& x );

private:
    int height_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalRow_;

    Matrix<T> localVec_;
};

// Set all of the entries of x to zero
template<typename T>
void MakeZeros( DistVector<T>& x );

// Draw the entries of x uniformly from the unitball in T
template<typename T>
void MakeUniform( DistVector<T>& x );

// Just an l2 norm for now
template<typename F>
typename Base<F>::type Norm( const DistVector<F>& x );

// y := alpha x + y
template<typename T>
void Axpy( T alpha, const DistVector<T>& x, DistVector<T>& y );

} // namespace cliq
