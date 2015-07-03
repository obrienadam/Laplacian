/**
*  
* This program solves the laplacian on a 2-D equidistant grid. It is meant
* to be an example of using petsc.
* 
* Author: Adam O'Brien
* Date  : 
*
**/

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <petsc.h>

typedef std::vector<double> FloatArray;
typedef std::vector<int> IntArray;

int main()
{
	using namespace std;

	int i, j, l, nI, nJ, nNodes, rowNo, nIters;
	double h, hSqr;
	double phiBE, phiBW, phiBN, phiBS, phi;
	string option, filename;
	FloatArray coeffs;
	IntArray jIndices;
	Mat A;
	Vec x, b;
	KSP ksp;
	PC pc;
	ofstream fout;

	cout << "Enter Number of Nodes in the X-direction : "; cin >> nI;
	cout << "Enter Number of Nodes in the Y-direction : "; cin >> nJ;
	cout << "Enter the Grid Node Spacing              : "; cin >> h;
	cout << "Enter east boundary condition value 	  : "; cin >> phiBE;
	cout << "Enter west boundary condition value 	  : "; cin >> phiBW;
	cout << "Enter north boundary condition value 	  : "; cin >> phiBN;
	cout << "Enter south boundary condition value 	  : "; cin >> phiBS;
	nNodes = nI*nJ;
	hSqr = h*h;
	phiBE /= -hSqr;
	phiBW /= -hSqr;
	phiBN /= -hSqr;
	phiBS /= -hSqr;

	coeffs.resize(5);
	jIndices.resize(5);

	PetscInitializeNoArguments();

	MatCreate(PETSC_COMM_SELF, &A);
	MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nNodes, nNodes);
	MatSetType(A, MATSEQAIJ);
	MatSeqAIJSetPreallocation(A, 5, PETSC_NULL);
	MatSetUp(A);

	VecCreate(PETSC_COMM_SELF, &b);
	VecSetSizes(b, PETSC_DECIDE, nNodes);
    VecSetType(b, VECSTANDARD);
    VecDuplicate(b, &x);
    VecZeroEntries(b);

	for(j = 0; j < nJ; ++j)
	{
		for(i = 0; i < nI; ++i)
		{
			coeffs[0] = -4./hSqr;
			coeffs[1] = 1./hSqr;
			coeffs[2] = 1./hSqr;
			coeffs[3] = 1./hSqr;
			coeffs[4] = 1./hSqr;

			rowNo = j*nI + i;
			jIndices[0] = rowNo;
			jIndices[1] = rowNo + 1;
			jIndices[2] = rowNo - 1;
			jIndices[3] = rowNo + nI;
			jIndices[4] = rowNo - nI;

			if(i == 0)
			{
				jIndices[2] = -1;
				VecSetValues(b, 1, &rowNo, &phiBW, ADD_VALUES);
			}
			if(i == nI - 1)
			{
				jIndices[1] = -1;
				VecSetValues(b, 1, &rowNo, &phiBE, ADD_VALUES);
			}
			if(j == 0)
			{
				jIndices[4] = -1;
				VecSetValues(b, 1, &rowNo, &phiBS, ADD_VALUES);
			}
			if(j == nJ - 1)
			{
				jIndices[3] = -1;
				VecSetValues(b, 1, &rowNo, &phiBN, ADD_VALUES);
			}

			MatSetValues(A, 1, &rowNo, 5, jIndices.data(), coeffs.data(), INSERT_VALUES);
		}
	}

	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	cout << "View matrix? (yes/no): "; cin >> option;
	std::transform(option.begin(), option.end(), option.begin(), ::tolower);

	if(option == "yes" || option == "y")
		MatView(A, PETSC_VIEWER_STDOUT_SELF);

	KSPCreate(PETSC_COMM_SELF, &ksp);
	KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCILU);
    PCFactorSetFill(pc, 2);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetType(ksp, KSPCG);
    KSPSolve(ksp, b, x);
    KSPGetIterationNumber(ksp, &nIters);

    cout << "The solution converged in " << nIters << " iterations." << endl
    	 << "Enter output filename : "; cin >> filename;

    fout.open(filename.c_str());

    for(j = 0; j < nJ; ++j)
    {
    	for(i = 0; i < nI; ++i)
    	{
    		rowNo = j*nI + i;
    		VecGetValues(x, 1, &rowNo, &phi);

    		fout << phi << ",";
    	}

    	fout << endl;
    }

	PetscFinalize();

	return 0;
}