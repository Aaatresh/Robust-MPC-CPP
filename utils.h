/*
    This file is a utility file in the robust MPC simulation and testing project. This file contains useful utility
    functions for use for mathematical operations.
*/

int sign(double x)
{
    /*
        This function implements the signum function. This function returns 1 if x is greater than 0.
        This function returns -1 if x is less than 0, and it returns 0 if x is equal to 0.

        Arguments:
            - x: Input number to signum.

        Returns:
            - {1, -1, 0} depending on the sign of x.
    */

    if(x > 0)
        return 1;
    else if(x < 0)
        return -1;
    else
        return 0;
}


double determinant(Matrix mat)
{
    /*
        This function determines the determinant of a matrix.

        This function was forked from https://www.geeksforgeeks.org/determinant-of-a-matrix/, and was modified by me
        to handle floating point numbers in the input matrix and to handle the matrix data structure that I am using.

        Arguments:
            - mat: Input matrix.

        Returns:
            - Determinant of the input matrix.
    */

    // Find the size of the input matrix
    int n = mat.i32getColumn();

    // Initialize result
    double num1, num2, det = 1, total = 1;
    int index;

    // Temporary array for storing row
    float temp[n + 1];

    // Loop for traversing the diagonal elements
    for (int i = 0; i < n; i++)
    {
        // Initialize the index
        index = i;

        // Finding the index which has
        // non zero value
        while (index < n && mat[index][i] == 0)
        {
            index++;
        }

        // if there is non zero element
        if (index == n)
        {
            // the determinant of matrix
            // as zero
            continue;
        }
        if (index != i)
        {
            // Loop for swapping the diagonal
            // element row and index row
            for (int j = 0; j < n; j++)
            {
                swap(mat[index][j], mat[i][j]);
            }

            // Determinant sign changes when we
            // shift rows go through determinant
            // properties
            det = det * pow(-1, index - i);
        }

        // Storing the values of diagonal
        // row elements
        for (int j = 0; j < n; j++)
        {
            temp[j] = mat[i][j];
        }

        // Traversing every row below the
        // diagonal element
        for (int j = i + 1; j < n; j++)
        {
            // Value of diagonal element
            num1 = temp[i];

            // Value of next row element
            num2 = mat[j][i];

            // Traversing every column of row
            // and multiplying to every row
            for (int k = 0; k < n; k++)
            {
                // Multiplying to make the diagonal
                // element and next row element equal
                mat[j][k]
                    = (num1 * mat[j][k]) - (num2 * temp[k]);
            }
            total = total * num1; // Det(kA)=kDet(A);
        }
    }

    // Multiplying the diagonal elements to
    // get determinant
    for (int i = 0; i < n; i++)
    {
        det = det * mat[i][i];
    }

    return (det / total);
}

double trace(Matrix M)
{
    /*
        Function to calculate the trace of an input matrix.

        Arguments:
            - M: Input matrix.

        Returns:
            - Trace of the matrix.
    */

    // Iterate through column indices and find the sum of all M[j][j].
    double sum = 0;
    for(int j=0; j < M.i32getColumn(); j++)
    {
        sum += M[j][j];
    }

    return sum;
}