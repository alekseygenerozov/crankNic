#ifndef INC_MR
#define INC_MR

/*
 *		M VECTOR
 *
 *			A small shell class to store an array
 *			with size information, smart intializers,
 *			assignment operator and optional bounds
 *			checking.
 */
template <class T> class Mvector {
private:
	size_t nn;
	T *v;
public:
	Mvector();                                  // constructors
	explicit Mvector( size_t n );
	Mvector(size_t n, const T &a );             // value
	Mvector(size_t n, const T *a );             // c-style array
	Mvector(const Mvector &rhs);                // copy
	Mvector & operator=(const Mvector &rhs);    // assignment
	typedef T val_type;
	inline T & operator[](const size_t i);              // indexing
	inline const T & operator[](const size_t i) const;  // "" for const
	inline size_t size() const;                 // figure out length
	void resize(size_t newn);                   // resize, losing contents
	void assign(size_t newn, const T &a);       // resize, assign a to spots
	~Mvector();
};

template <class T> Mvector<T>::Mvector()
 : nn(0), v(NULL) {}

template <class T> Mvector<T>::Mvector(size_t n)
 : nn(n), v( n > 0 ? new T[n] : NULL ) {}

template <class T> Mvector<T>::Mvector(size_t n, const T &a)
 : nn(n), v( n > 0 ? new T[n] : NULL ) 
{
	for( T *i = v ; i != v + nn ; ++i )
		*i = a;
}

template <class T> Mvector<T>::Mvector( size_t n, const T *a )
 : nn(n) , v( n > 0 ? new T[n] : NULL )
{
	for( T *lhs = v, *rhs = a; lhs != v + nn ; ++lhs, ++rhs )
			*lhs = *rhs;
}

template <class T> Mvector<T>::Mvector(const Mvector<T> &rhs)
 : nn(rhs.nn) , v( nn > 0 ? new T[nn] : NULL )
{
	for( size_t i = 0 ; i < nn ; ++i )
		v[i] = rhs[i];
}

template <class T> Mvector<T> &Mvector<T>::operator=(const Mvector<T> &rhs )
{
	if( this != &rhs )
	{
		if( nn != rhs.nn ) {
			if( v != NULL ) delete [] (v);
			nn = rhs.nn;
			v = nn > 0 ? new T[nn] : NULL;
		}
		for( size_t i = 0 ; i < nn ; ++i )
			v[i] = rhs[i];
	}
	return *this;
}

template <class T> inline T &Mvector<T>::operator[](const size_t i)
{
#ifdef _CHECKBNDS_
	if( i < 0 || i >= nn ) throw("vector index out of bounds");
#endif
	return v[i];
}

template <class T> inline const T &Mvector<T>::operator[](const size_t i) const
{
#ifdef _CHECKBNDS_
	if( i < 0 || i >= nn ) throw("vector index out of bounds");
#endif
	return v[i];
}

template <class T> inline size_t Mvector<T>::size() const 
{ return nn; }

template <class T> void Mvector<T>::resize(size_t newn)
{
	if( newn != nn ){
		if( v != NULL) delete [] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T> void Mvector<T>::assign(size_t newn, const T &a)
{
	if( newn != nn ){
		if( v != NULL ) delete [] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for( T *i = v ; i != v + nn ; ++i )
		*i = a;
}

template <class T> Mvector<T>::~Mvector()
{
	if( v != NULL ) delete [] (v);
}


/*
 *	M MATRIX
 *
 *
 *
 */
template <class T> class Mmatrix {
private:
	size_t nn;
	size_t mm;
	T **v;
public:
	Mmatrix();
	Mmatrix(size_t n, size_t m);
	Mmatrix(size_t n, size_t m, const T &a);
	Mmatrix(size_t n, size_t m, const T *a);
	Mmatrix(const Mmatrix &rhs);
	Mmatrix & operator=( const Mmatrix &rhs );
	typedef T val_type;
	inline T* operator[](const size_t i);
	inline const T* operator[]( const size_t i ) const;
	inline size_t nrows() const;
	inline size_t ncols() const;
	void resize(size_t n, size_t m);
	void assign( size_t n, size_t m, const T &a);
	~Mmatrix();
};// end Mmatrix

template <class T> Mmatrix<T>::Mmatrix() : nn(0), mm(0), v(NULL) {}

template <class T> Mmatrix<T>::Mmatrix( size_t n , size_t m )
	: nn(n) , mm(m) , v( n > 0 ? new T*[n] : NULL )
{
	size_t i, nbm = m*n;
	if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
	for( i = 1 ; i < n ; ++i ) v[i] = v[i-1] + m;
} // end constructor

template <class T> Mmatrix<T>::Mmatrix( size_t n , size_t m , const T &a )
	: nn(n) , mm(m) , v( n > 0 ? new T*[n] : NULL )
{
	size_t i,j,nbm=mm*nn;
	if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
	for( i = 1 ; i < n ; ++i ) v[i] = v[i-1] + m;
	for( i = 0 ; i < n ; ++i ) for( j = 0 ; j < m ; ++j ) v[i][j] = a;
} // end constructor

template <class T> Mmatrix<T>::Mmatrix( size_t n , size_t m , const T *a )
	: nn(n) , mm(m) , v( n > 0 ? new T*[n] : NULL )
{
	size_t i,j,nbm=m*n;
	if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
	for( i = 1 ; i < n ; ++i ) v[i] = v[i-1] + m;
	for( i = 0 ; i < n ; ++i ) for( j = 0 ; j < m ; ++j ) v[i][j] = *a++;
} // end constructor

template <class T> Mmatrix<T>::Mmatrix( const Mmatrix<T> &rhs )
	: nn(rhs.nn), mm(rhs.mm) , v( nn > 0 ? new T*[nn] : NULL )
{
	size_t i,j,nbm=mm*nn;
	if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
	for( i = 1 ; i < nn ; ++i ) v[i] = v[i-1] + mm;
	for( i = 0 ; i < nn ; ++i ) for( j = 0 ; j < mm ; ++j ) v[i][j] = rhs[i][j];
}// end constructor

template <class T> Mmatrix<T> & Mmatrix<T>::operator=(const Mmatrix<T> &rhs)
{
	if( this != &rhs ){
		size_t i,j,nbm;
		if( nn != rhs.nn || mm != rhs.mm ){
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn > 0 ? new T*[nn] : NULL;
			nbm = mm*nn;
			if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
			for( i = 1 ; i < nn ; ++i ) v[i] = v[i-1] + mm;
		}// if dimensions don't match
		for( i = 0 ; i < nn ; ++i ) for ( j = 0 ; j < mm ; j++ ) v[i][j] = rhs[i][j];
	}// end check if
	return *this;
}// end = op

template <class T> inline T* Mmatrix<T>::operator[]( const size_t i )
{
#ifdef _CHECKBNDS_
if( i < 0 || i >= nn ) throw("Matrix subscript out of bounds!");
#endif
	return v[i];
}// end [] op

template <class T> inline const T* Mmatrix<T>::operator[]( const size_t i ) const {
#ifdef _CHECKBNDS_
if(i<0 || i >= nn) throw("Matrix subscript out of bounds!");
#endif
	return v[i];
}// end [] op

template <class T> inline size_t Mmatrix<T>::nrows() const { return nn; }
template <class T> inline size_t Mmatrix<T>::ncols() const { return mm; }

template <class T> void Mmatrix<T>::resize( size_t n, size_t m )
{
	size_t i,nbm;
	if( n != nn || m != mm ){
		if( v != NULL ){
			delete [] (v[0]);
			delete [] (v);
		}// end NULL if
		nn = n;
		mm = m;
		v = nn > 0 ? new T*[n] : NULL;
		nbm = nn*mm;
		if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
		for( i = 1 ; i < nn ; ++i ) v[i] = v[i-1] + mm;
	}// end dims if
}// end resize

template <class T> void Mmatrix<T>::assign( size_t n , size_t m , const T& a)
{
	size_t i,j,nbm;
	if( n != nn || m != mm ){
		if( v != NULL ) {
			delete [] (v[0]);
			delete [] (v);
		} // end NULL if
		nn = n;
		mm = m;
		v = nn > 0 ? new T*[nn] : NULL;
		nbm = mm*nn;
		if(v) v[0] = nbm > 0 ? new T[nbm] : NULL;
		for( i = 1 ; i < nn ; ++i ) v[i] = v[i-1] + mm;
	}// end dim if
	for( i = 0 ; i < nn ; ++i ) for( j = 0 ; j < mm ; ++j ) v[i][j] = a;
}// end assign

template <class T> Mmatrix<T>::~Mmatrix()
{
	if( v != NULL ){
		delete [] (v[0]);
		delete [] (v);
	}// end NULL if
}// end destructor



// Vector typedefs

typedef const Mvector<int> vInt_i;
typedef Mvector<int> vInt, vInt_o, vInt_io;

typedef const Mvector<double> vDoub_i;
typedef Mvector<double> vDoub, vDoub_o, vDoub_io;

typedef const Mvector<double*> vDoubp_i;
typedef Mvector<double*> vDoubp, vDoubp_o, vDoubp_io;

// Matrix typedefs

typedef const Mmatrix<int> mInt_i;
typedef Mmatrix<int> mInt, mInt_o, mInt_io;

typedef const Mmatrix<double> mDoub_i;
typedef Mmatrix<double> mDoub, mDoub_o, mDoub_io;

#ifndef _NR3_H_			// FIXME
// Useful helper function
template<class T> inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
#endif

#endif
