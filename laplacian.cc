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
#include <mpi.h>
#include <petsc.h>

int main(int argc, char *argv[])
{
  using namespace std;

  int nI = 100, nJ = 100;
  double h = 0.01, phiE = 0, phiW = 1, phiN = 0, phiS = 0;

  MPI_Init(&argc, &argv);
  PetscInitializeNoArguments();

  for(int i = 0; i < argc; ++i)
    {
      string arg(argv[i]);

      if(arg == "--n")
	nI = nJ = stoi(argv[++i]);
      else if(arg == "--h")
	h = stod(argv[++i]);
      else if(arg == "--phiE")
	phiE = stod(argv[++i]);
      else if(arg == "--phiW")
	phiW = stod(argv[++i]);
      else if(arg == "--phiN")
	phiN = stod(argv[++i]);
      else if(arg == "--phiS")
	phiS = stod(argv[++i]);
    }

  const int nNodes = nI*nJ;
  const double hSqr = h*h;
  const double coeffs[] = {1/hSqr, 1/hSqr, -4/hSqr, 1/hSqr, 1/hSqr};

  phiE /= -hSqr;
  phiW /= -hSqr;
  phiN /= -hSqr;
  phiS /= -hSqr;

  Mat A;
  Vec x, b;
  MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nNodes, nNodes, 5, NULL, 3, NULL, &A);
  pair<int, int> range;
  MatGetOwnershipRange(A, &range.first, &range.second);

  

  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nNodes, &b);
  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, nNodes, &x);
  VecZeroEntries(b);

  for(int j = 0; j < nJ; ++j)
    for(int i = 0; i < nI; ++i)
      {
	int k = nI*j + i;

	if(k < range.first)
	  continue;
	else if(k >= range.second)
	  goto endloop;

	int cols[] = {k - nI, k - 1, k, k + 1, k + nI};

	if(i == 0)
	 {
	   VecSetValues(b, 1, &k, &phiW, ADD_VALUES);
	   cols[1] = -1;
	 }
	if(i == nI - 1)
	  {
	    VecSetValues(b, 1, &k, &phiE, ADD_VALUES);
	    cols[3] = -1;
	  }
	if(j == 0)
	  {
	    VecSetValues(b, 1, &k, &phiS, ADD_VALUES);
	    cols[0] = -1;
	  }
	if(j == nJ - 1)
	  {
	    VecSetValues(b, 1, &k, &phiN, ADD_VALUES);
	    cols[4] = -1;
	  }
	
	MatSetValues(A, 1, &k, 5, cols, coeffs, INSERT_VALUES);
      }
  endloop:

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  KSP ksp;
  PC pc;
  int nIters;

  KSPCreate(MPI_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCSOR);
  KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetType(ksp, KSPCG);
  KSPSolve(ksp, b, x);
  KSPGetIterationNumber(ksp, &nIters);

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0)
    {
      ofstream fout;
      fout.open("result.txt");

    for(int j = 0; j < nJ; ++j)
      {
        for(int i = 0; i < nI; ++i)
    	  {
	    int k = j*nI + i;
	    if(k < range.first)
	      continue;
	    else if(k >= range.second)
	      goto endloop2;

	    double phi;
	    VecGetValues(x, 1, &k, &phi);

	    fout << phi << ",";
    	  }
      
        fout << endl;
      }
      endloop2:
      fout.close();
    }

  MatDestroy(&A);
  PetscFinalize();
  MPI_Finalize();

  return 0;
}
