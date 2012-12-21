/*
   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This file is part of Clique and is under the GNU General Public License,
   which can be found in the LICENSE file in the root directory, or at 
   <http://www.gnu.org/licenses/>.
*/

namespace cliq {

// Use a simple 1d distribution where each process owns a fixed number of 
// indices,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
class DistMap
{
public:
    // Constructors and destructors
    DistMap();
    DistMap( mpi::Comm comm );
    DistMap( int numSources, mpi::Comm comm );
    // TODO: Constructor for building from a DistMap
    ~DistMap();

    // If the total number of sources is partitioned among the processes, 
    // calling this routine will have the DistMap map the indices to the 
    // owning process
    void StoreOwners
    ( int numSource, std::vector<int>& localIndices, mpi::Comm comm );

    // Map manipulation
    // Collectively map each process's local set of indices
    void Translate( std::vector<int>& localIndices ) const;
    // Form the inverse map
    void FormInverse( DistMap& inverseMap ) const;
    // composite(i) := second(first(i))
    void Extend( DistMap& firstMap ) const;
    void Extend( const DistMap& firstMap, DistMap& compositeMap ) const;

    // High-level information
    int NumSources() const;

    // Communicator management
    void SetComm( mpi::Comm comm );
    mpi::Comm Comm() const;

    // Distribution information
    int Blocksize() const;
    int FirstLocalSource() const;
    int NumLocalSources() const;

    // Local data
    int GetLocal( int localSource ) const;
    void SetLocal( int localSource, int target );
    int* LocalBuffer();
    const int* LocalBuffer() const;
    std::vector<int>& LocalMap();
    const std::vector<int>& LocalMap() const;

    // For modifying the size of the map
    void Empty();
    void ResizeTo( int numSources );

    // Assignment
    const DistMap& operator=( const DistMap& map );

private:
    int numSources_;

    mpi::Comm comm_;

    int blocksize_;
    int firstLocalSource_;

    std::vector<int> localMap_;
};

} // namespace cliq
