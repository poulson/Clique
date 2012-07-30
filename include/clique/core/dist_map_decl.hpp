/*
   Clique: a scalable implementation of the multifrontal algorithm

   Copyright (C) 2011-2012 Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
