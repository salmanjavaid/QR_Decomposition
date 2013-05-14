#include "QR_Header.hpp"


/* QR decomposition using Householder Reflectors wiki link: http://en.wikipedia.org/wiki/Householder_transformation */


/* Helper function to convert a matrix to a vector */

vector<double> ConvertMatrix2Vector(const matrix<double> &input)
{
	
	vector<double> R(input.size1() * input.size2());    /* create a new vector of size m*n  of input */
	int Index = 0;
	for(int i = 0; i < input.size1(); i++)              /* run a loop and update each row to the newly created vector */
	{
		for(int j = 0; j < input.size2(); j++)
		{
			R(Index) =  input(i,j);
			Index++;
		}
	}
	return R;
}


void QR(matrix<double> &input, matrix<double> &Q, matrix<double> &R)
{
	size_t Rows = input.size1();    /* Row size of input matrix */
	size_t Cols = input.size2();    /* Col size of input matrix */

	identity_matrix<double> I (Rows, Rows);  /* create an identity matrix for Householder Matrix computation */

	/* helper matrices and vectors */

	matrix<double> temp_Mat;
	vector<double> Vec;
	vector<double> d;
	

	for(int i = 0; i < (int)(std::min(Rows - 1, Cols)); i++)            
	{
		
		matrix<double> H;    /* helper matrix to store Householder matrix */
		Vec = ConvertMatrix2Vector(subrange(input, i, Rows, i, i + 1));  /* Take the column according to the iteration and convert to vector */
		temp_Mat = subrange(input, i, Rows, i, Cols); /* helper matrix */
		
		if (temp_Mat(0,0) > 0)     /* depending upon the value at first index decide to subtract or add vector to compute reflector*/
		{
			Vec = Vec + (norm_2(Vec) * ConvertMatrix2Vector(subrange(I, i, Rows, i, i + 1)));
		}
		else
		{
			Vec = Vec - (norm_2(Vec) * ConvertMatrix2Vector(subrange(I, i, Rows, i, i + 1)));
		}

		

		matrix<double> V(Vec.size(), 1);   /* uBlas supports matrix to matrix multiplication thus convert Vector to matrix */
		V <<= Vec;  /* assign vector to matrix */
		double n_1 = norm_1(prod(trans(V), V));  /* this quantity will always be a scalar but prod returns a vector so take norm */
												 /* to convert it to double */

		if (n_1 > 0.00025)           /* this value is often very small thus don't evaluate following code if its too small. 0.00025 is heuristic choice */
		{

			H = subrange(I, i, Rows, i, Rows) - (double)(2) * prod(V, trans(V)) / n_1;  /* calculate householder matrix */
			subrange(input, i, Rows, i, Cols) = prod(H, subrange(input, i, Rows, i, Cols)); /* calculate A(i) */
		
		
			/* Steps to compute Q, in this case we multiply our original matrix of size m*m by m - i * m - i  matrices successively */
			/* the rest of m - i * m - i matrix will be filled with Identity matrix */

			if (i > 0)				/* if i is greater than zero then we need a helper identity matrix */
			{
				identity_matrix<double> T_1 (Rows, Rows);     /* create helper identity matrix */
				matrix<double> _P = T_1;					  /* copy the matrix into a helper matrix */
				subrange(_P, i, Rows, i, Rows) = H;			  /* copy householder matrix into helper matrix */
				Q = prod(Q, _P);							  /* Calculate Q at this iteration */
			}
			else
			{
				Q = H;           /* at the first step no need of helper matrix, just copy Q into H */
			}
		}
	}

	R = input;   /* copy input to R */
}



int main () 
{
    
    matrix<double> m (4, 2);  /* this is our input matrix  */
	matrix<double> q (m.size1(), m.size1());  /* this is our Q matrix */
	matrix<double> r (m.size1(), m.size2());  /* this is our R matrix */


	/* populate m */

	int Index = 1;
    for (int i = 0; i < (int)(m.size1 ()); ++ i)
	{
		
        for (int j = 0; j < (int)(m.size2 ()); ++ j)
		{
            m (i, j) = (double)(Index);
			Index++;
		}
	}
    
	/*clock_t start, end;
	double msecs;
	start = clock();*/
	
	
	QR(m, q, r);   /* calling QR for decomposition */
	std::cout << q << std::endl;  /* print q */
	std::cout << r << std::endl;  /* print r */

	
	/*

	end = clock();
	msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
	std::cout << msecs << std::endl;


	*/

	return 0;
}