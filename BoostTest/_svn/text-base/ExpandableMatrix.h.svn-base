#ifndef EXPANDABLE_MATRIX
#define EXPANDABLE_MATRIX

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include	<boost/function.hpp>

using namespace boost::numeric::ublas;


class ExpandableMatrixPair
{
	matrix<double> _horz;
	matrix<double> _vert;

	boost::numeric::ublas::range _horzrange, _vertrange;

	int _vertColCount,		//Number of columns in the vertical matrix
		_horzRowCount,		//Number of rows in the horizontal matrix
		_size,				//statecount, #cols of horz = #rows of vert
		_totalSize;			// _matrixCount * _size

	int _horzIndex, _vertIndex;

	double _dt;

public:
	ExpandableMatrixPair(int matrixCount, int size, int horzRowCount, int vertColCount, double dt);
	void Add(matrix< double>& horz, matrix< double>& vert, matrix<double>& temp);


};
#endif