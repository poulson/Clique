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
class Vector
{
public:
    // Constructors and destructors
    Vector();
    Vector( int height );
    // TODO: Constructor for building from a Vector
    ~Vector();

    // High-level information
    int Height() const;

    // Data
    T Get( int row ) const;
    void Set( int row, T value );
    void Update( int row, T value );

    // For modifying the size of the vector
    void Empty();
    void ResizeTo( int height );

    // Assignment
    const Vector<T>& operator=( const Vector<T>& x );

private:
    Matrix<T> vec_;
};

// Set all of the entries of x to zero
template<typename T>
void MakeZeros( Vector<T>& x );

// Draw the entries of x uniformly from the unitball in T
template<typename T>
void MakeUniform( Vector<T>& x );

// Just an l2 norm for now
template<typename F>
typename Base<F>::type Norm( const Vector<F>& x );

// y := alpha x + y
template<typename T>
void Axpy( T alpha, const Vector<T>& x, Vector<T>& y );

} // namespace cliq
