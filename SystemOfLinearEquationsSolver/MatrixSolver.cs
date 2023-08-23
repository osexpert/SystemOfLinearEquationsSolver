﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace SystemOfLinearEquationsSolver
{
	/// <summary>
	/// https://stackoverflow.com/a/46837817/2671330
	/// https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
	/// 
	/// </summary>
	public class MatrixSolver_LU : ISolver
	{
		public static readonly MatrixSolver_LU Instance = new MatrixSolver_LU();

		private MatrixSolver_LU()
		{

		}

		public double[] SolveEquations(double[,] coefficients, int? resultRoundDigits = null)
		{
			var dim0Length = coefficients.GetLength(0);
			var dim1Length = coefficients.GetLength(1);

			//if (abCombined.Length != abCombined[0].Length - 1)
			if (dim0Length != dim1Length - 1)
				throw new SolverException("Last should be the b-matrix");

			double[][] m = MatrixCreate(dim0Length, dim0Length);
			double[] b = new double[dim0Length];
			for (int i = 0; i < dim0Length; i++)
			{
				//				double[] row = abCombined[i];
				//for (int j = 0; j < row.Length - 1; j++)
				for (int j = 0; j < dim1Length - 1; j++)
				{
					m[i][j] = coefficients[i, j];
				}

				//b[i] = -row[row.Length - 1]; //PS: negative!
				b[i] = coefficients[i, dim1Length - 1];
			}

			var resDec = SystemSolve(m, b);

			if (resultRoundDigits != null)
			{
				for (int k = 0; k < resDec.Length; k++)
				{
					resDec[k] = Math.Round(resDec[k], resultRoundDigits.Value);
				}
			}

			return resDec;
		}

		static double[] SystemSolve(double[][] A, double[] b)
		{
			// Solve Ax = b
			int n = A.Length;

			// 1. decompose A
			int[] perm;
			int toggle;
			double[][] luMatrix = MatrixDecompose(A, out perm, out toggle);
			if (luMatrix == null)
				return null;

			// 2. permute b according to perm[] into bp
			double[] bp = new double[b.Length];
			for (int i = 0; i < n; ++i)
				bp[i] = b[perm[i]];

			// 3. call helper
			double[] x = HelperSolve(luMatrix, bp);
			return x;
		} // SystemSolve

		static double[][] MatrixCreate(int rows, int cols)
		{
			double[][] result = new double[rows][];
			for (int i = 0; i < rows; ++i)
				result[i] = new double[cols];
			return result;
		}

		static double[][] MatrixIdentity(int n)
		{
			// return an n x n Identity matrix
			double[][] result = MatrixCreate(n, n);
			for (int i = 0; i < n; ++i)
				result[i][i] = 1d;

			return result;
		}

		static double[][] MatrixProduct(double[][] matrixA, double[][] matrixB)
		{
			int aRows = matrixA.Length; int aCols = matrixA[0].Length;
			int bRows = matrixB.Length; int bCols = matrixB[0].Length;
			if (aCols != bRows)
				throw new SolverException("Non-conformable matrices in MatrixProduct");

			double[][] result = MatrixCreate(aRows, bCols);

			for (int i = 0; i < aRows; ++i) // each row of A
				for (int j = 0; j < bCols; ++j) // each col of B
					for (int k = 0; k < aCols; ++k) // could use k less-than bRows
						result[i][j] += matrixA[i][k] * matrixB[k][j];

			return result;
		}

		static double[][] MatrixInverse(double[][] matrix)
		{
			int n = matrix.Length;
			double[][] result = MatrixDuplicate(matrix);

			int[] perm;
			int toggle;
			double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
			if (lum == null)
				throw new Exception("Unable to compute inverse");

			double[] b = new double[n];
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					if (i == perm[j])
						b[j] = 1d;
					else
						b[j] = 0d;
				}

				double[] x = HelperSolve(lum, b);

				for (int j = 0; j < n; ++j)
					result[j][i] = x[j];
			}
			return result;
		}

		static double[][] MatrixDuplicate(double[][] matrix)
		{
			// allocates/creates a duplicate of a matrix.
			double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
			for (int i = 0; i < matrix.Length; ++i) // copy the values
				for (int j = 0; j < matrix[i].Length; ++j)
					result[i][j] = matrix[i][j];
			return result;
		}

		static double[] HelperSolve(double[][] luMatrix, double[] b)
		{
			// before calling this helper, permute b using the perm array
			// from MatrixDecompose that generated luMatrix
			int n = luMatrix.Length;
			double[] x = new double[n];
			b.CopyTo(x, 0);

			for (int i = 1; i < n; ++i)
			{
				double sum = x[i];
				for (int j = 0; j < i; ++j)
					sum -= luMatrix[i][j] * x[j];
				x[i] = sum;
			}

			x[n - 1] /= luMatrix[n - 1][n - 1];
			for (int i = n - 2; i >= 0; --i)
			{
				double sum = x[i];
				for (int j = i + 1; j < n; ++j)
					sum -= luMatrix[i][j] * x[j];
				x[i] = sum / luMatrix[i][i];
			}

			return x;
		}

		static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
		{
			// Doolittle LUP decomposition with partial pivoting.
			// rerturns: result is L (with 1s on diagonal) and U;
			// perm holds row permutations; toggle is +1 or -1 (even or odd)
			int rows = matrix.Length;
			int cols = matrix[0].Length; // assume square
			if (rows != cols)
				throw new Exception("Attempt to decompose a non-square m");

			int n = rows; // convenience

			double[][] result = MatrixDuplicate(matrix);

			perm = new int[n]; // set up row permutation result
			for (int i = 0; i < n; ++i) { perm[i] = i; }

			toggle = 1; // toggle tracks row swaps.
						// +1 -greater-than even, -1 -greater-than odd. used by MatrixDeterminant

			for (int j = 0; j < n - 1; ++j) // each column
			{
				double colMax = Math.Abs(result[j][j]); // find largest val in col
				int pRow = j;
				//for (int i = j + 1; i less-than n; ++i)
				//{
				//  if (result[i][j] greater-than colMax)
				//  {
				//    colMax = result[i][j];
				//    pRow = i;
				//  }
				//}

				// reader Matt V needed this:
				for (int i = j + 1; i < n; ++i)
				{
					var abs = Math.Abs(result[i][j]);
					if (abs > colMax)
					{
						colMax = abs;// Math.Abs(result[i][j]);
						pRow = i;
					}
				}
				// Not sure if this approach is needed always, or not.

				if (pRow != j) // if largest value not on pivot, swap rows
				{
					double[] rowPtr = result[pRow];
					result[pRow] = result[j];
					result[j] = rowPtr;

					int tmp = perm[pRow]; // and swap perm info
					perm[pRow] = perm[j];
					perm[j] = tmp;

					toggle = -toggle; // adjust the row-swap toggle
				}

				// --------------------------------------------------
				// This part added later (not in original)
				// and replaces the 'return null' below.
				// if there is a 0 on the diagonal, find a good row
				// from i = j+1 down that doesn't have
				// a 0 in column j, and swap that good row with row j
				// --------------------------------------------------

				if (result[j][j] == 0d)
				{
					// find a good row to swap
					int goodRow = -1;
					for (int row = j + 1; row < n; ++row)
					{
						if (result[row][j] != 0d)
							goodRow = row;
					}

					if (goodRow == -1)
						throw new SolverException("Cannot use Doolittle's method");

					// swap rows so 0.0 no longer on diagonal
					double[] rowPtr = result[goodRow];
					result[goodRow] = result[j];
					result[j] = rowPtr;

					int tmp = perm[goodRow]; // and swap perm info
					perm[goodRow] = perm[j];
					perm[j] = tmp;

					toggle = -toggle; // adjust the row-swap toggle
				}
				// --------------------------------------------------
				// if diagonal after swap is zero . .
				//if (Math.Abs(result[j][j]) less-than 1.0E-20) 
				//  return null; // consider a throw

				for (int i = j + 1; i < n; ++i)
				{
					result[i][j] /= result[j][j];
					for (int k = j + 1; k < n; ++k)
					{
						result[i][k] -= result[i][j] * result[j][k];
					}
				}


			} // main j column loop

			return result;
		}

	
	}
}
