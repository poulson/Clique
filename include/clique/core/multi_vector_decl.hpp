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

template<typename T>
class MultiVector
{
public:
    // Constructors and destructors
    MultiVector();
    MultiVector( int height, int width );
    // TODO: Constructor for building from a MultiVector
    ~MultiVector();

    // High-level information
    int Height() const;
    int Width() const;

    // Data
    T Get( int row, int col ) const;
    void Set( int row, int col, T value );
    void Update( int row, int col, T value );

    // For modifying the size of the multi-vector
    void Empty();
    void ResizeTo( int height, int width );

    // Assignment
    const MultiVector<T>& operator=( const Vector<T>& x );
    const MultiVector<T>& operator=( const MultiVector<T>& X );

private:
    Matrix<T> multiVec_;
};

// Set all of the entries of X to zero
template<typename T>
void MakeZeros( MultiVector<T>& X );

// Draw the entries of X uniformly from the unitball in T
template<typename T>
void MakeUniform( MultiVector<T>& X );

// Just an l2 norm for now
template<typename F>
void Norms
( const MultiVector<F>& X, std::vector<typename Base<F>::type>& norms );

// Y := alpha X + Y 
template<typename T>
void Axpy( T alpha, const MultiVector<T>& X, MultiVector<T>& Y );

} // namespace cliq
