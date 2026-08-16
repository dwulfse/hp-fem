#include "CSR_Matrix.hpp"
#include <stdexcept>
#include <iostream>

// constructors
CSRMatrix::CSRMatrix()
{
}

CSRMatrix::CSRMatrix(int nnz, int n)
 : entries(nnz), col_no(nnz), row_start(n+1)
{
}

CSRMatrix::CSRMatrix(std::vector<double> entries, std::vector<int> col_no, std::vector<int> row_start)
 : entries(entries), col_no(col_no), row_start(row_start)
{
}

// destructor
CSRMatrix::~CSRMatrix()
{
}

// override () operator for read-only indexing
double& CSRMatrix::operator()(int i, int j)
{
	int noRows = static_cast<int>(row_start.size()) - 1;
	if (i >= noRows || i < 0 || j >= noRows || j < 0)
	{
		throw std::out_of_range("Index out of range");
	}

	for (int k=row_start[i]; k<row_start[i+1]; k++)
	{
		if (col_no[k] == j)
		{
			return entries[k];
		}
	}

	throw std::out_of_range("Element not accessible");
}

void CSRMatrix::print(bool sparcity)
{
    int noRows = row_start.size() - 1;

    if (sparcity)
    {
        for (int i = 0; i < noRows; i++)
        {
            int k = row_start[i];
            for (int j = 0; j < noRows; j++)
            {
                if (k < row_start[i+1] && col_no[k] == j)
                {
                    std::cout << "X ";
                    k++;
                }
                else
                {
                    std::cout << "- ";
                }
            }
            std::cout << "\n";
        }
    }
    else
    {
        for (int i = 0; i < noRows; i++)
        {
            int k = row_start[i];
            for (int j = 0; j < noRows; j++)
            {
                if (k < row_start[i+1] && col_no[k] == j)
                {
                    std::cout << entries[k] << "\t ";
                    k++;
                }
                else
                {
                    std::cout << 0 << "\t ";
                }
            }
            std::cout << "\n";
        }
    }
}