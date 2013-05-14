
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <ctime>
using namespace boost::numeric::ublas;


void QR(matrix<double> &input, matrix<double> &Q, matrix<double> &R);
vector<double> ConvertMatrix2Vector(matrix<double> &input);