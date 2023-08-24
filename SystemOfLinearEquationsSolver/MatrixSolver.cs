using System;
using System.Collections.Generic;
using System.Linq;

namespace SystemOfLinearEquationsSolver
{
	/// <summary>
	/// https://stackoverflow.com/a/46837817/2671330
	/// https://jamesmccaffrey.wordpress.com/2015/03/06/inverting-a-matrix-using-c/
	/// 
	/// 
	///  http://msdn.microsoft.com/en-us/magazine/jj863137.aspx
	/// library of static matrix methods including decomposition, inverse, determinant
	/// Anonymous9807: Friday, February 1, 2013 3:48 AM
	/// "There appears to be a bug in the MatrixDecompose method. On lines 17 and 19, I think you should use Math.Abs() around the matrix terms "result[i][j]". 
	/// I suspect you will run into trouble when a matrix row has 0 on the main diagonal and negative values otherwise."
	/// 
	/// 
	/// </summary>
	public class MatrixSolver : ISolver
	{
		public static readonly MatrixSolver Instance = new MatrixSolver();

		private MatrixSolver()
		{

		}

		public double[] SolveEquations(double[,] coefficients, int? resultRoundDigits = null)
		{
			var numRows = coefficients.GetLength(0);
			var numCols = coefficients.GetLength(1);

			if (numRows != numCols - 1)
				throw new SolverException("Last should be the b-matrix (assume square A-matrix | B-matrix)");

			double[][] m = MatrixCreate(numRows, numRows);
			double[] b = new double[numRows];
			for (int i = 0; i < numRows; i++)
			{
				//				double[] row = abCombined[i];
				//for (int j = 0; j < row.Length - 1; j++)
				for (int j = 0; j < numCols - 1; j++)
				{
					m[i][j] = coefficients[i, j];
				}

				b[i] = coefficients[i, numCols - 1];
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
			// allocates/creates a matrix initialized to all 0.0. assume rows and cols > 0
			// do error checking here
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
					for (int k = 0; k < aCols; ++k) // could use k < bRows
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
			//double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
			//for (int i = 0; i < matrix.Length; ++i) // copy the values
			//	for (int j = 0; j < matrix[i].Length; ++j)
			//		result[i][j] = matrix[i][j];
			//return result;
			return (double[][])matrix.Clone();
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

		/// <summary>
		/// This became twize as slow when rewrite to jagged arrays, because row swap became very expensive.
		/// Even BlockCopy did not help.
		/// </summary>
		static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
		{
			// Doolittle LUP decomposition with partial pivoting.
			// rerturns: result is L (with 1s on diagonal) and U; perm holds row permutations; toggle is +1 or -1 (even or odd)
			int rows = matrix.Length;
			int cols = matrix[0].Length; // assume square (assume all rows have the same number of columns so just use row [0])
			if (rows != cols)
				throw new Exception("Attempt to decompose a non-square matrix");

			int n = rows; // convenience

			double[][] result = MatrixDuplicate(matrix);

			perm = new int[n]; // set up row permutation result
			for (int i = 0; i < n; ++i) { perm[i] = i; }

			toggle = 1; // toggle tracks row swaps. +1 -> even, -1 -> odd. used by MatrixDeterminant

			for (int j = 0; j < n - 1; ++j) // each column
			{
				double colMax = Math.Abs(result[j][j]); // find largest val in col j
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
					// Swap rows
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


		static double[][] MatrixRandom(int rows, int cols, double minVal, double maxVal, int seed)
		{
			// return a matrix with random values
			Random ran = new(seed);
			double[][] result = MatrixCreate(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					result[i][j] = (maxVal - minVal) * ran.NextDouble() + minVal;
			return result;
		}


		static string MatrixAsString(double[][] matrix)
		{
			string s = "";
			for (int i = 0; i < matrix.Length; ++i)
			{
				for (int j = 0; j < matrix[i].Length; ++j)
					s += matrix[i][j].ToString("F3").PadLeft(8) + " ";
				s += Environment.NewLine;
			}
			return s;
		}

		static bool MatrixAreEqual(double[][] matrixA, double[][] matrixB, double epsilon)
		{
			// true if all values in matrixA == corresponding values in matrixB
			int aRows = matrixA.Length; int aCols = matrixA[0].Length;
			int bRows = matrixB.Length; int bCols = matrixB[0].Length;
			if (aRows != bRows || aCols != bCols)
				throw new Exception("Non-conformable matrices in MatrixAreEqual");

			for (int i = 0; i < aRows; ++i) // each row of A and B
				for (int j = 0; j < aCols; ++j) // each col of A and B												
					if (Math.Abs(matrixA[i][j] - matrixB[i][j]) > epsilon)
						return false;
			return true;
		}


		static double[] MatrixVectorProduct(double[][] matrix, double[] vector)
		{
			// result of multiplying an n x m matrix by a m x 1 column vector (yielding an n x 1 column vector)
			int mRows = matrix.Length; int mCols = matrix[0].Length;
			int vRows = vector.Length;
			if (mCols != vRows)
				throw new Exception("Non-conformable matrix and vector in MatrixVectorProduct");
			double[] result = new double[mRows]; // an n x m matrix times a m x 1 column vector is a n x 1 column vector
			for (int i = 0; i < mRows; ++i)
				for (int j = 0; j < mCols; ++j)
					result[i] += matrix[i][j] * vector[j];
			return result;
		}

		static double MatrixDeterminant(double[][] matrix)
		{
			double[][] lum = MatrixDecompose(matrix, out int[] perm, out int toggle);
			if (lum == null)
				throw new Exception("Unable to compute MatrixDeterminant");
			double result = toggle;
			for (int i = 0; i < lum.Length; ++i)
				result *= lum[i][i];
			return result;
		}

		static double[][] ExtractLower(double[][] matrix)
		{
			// lower part of a Doolittle decomposition (1.0s on diagonal, 0.0s in upper)
			int rows = matrix.Length; int cols = matrix[0].Length;
			double[][] result = MatrixCreate(rows, cols);
			for (int i = 0; i < rows; ++i)
			{
				for (int j = 0; j < cols; ++j)
				{
					if (i == j)
						result[i][j] = 1.0;
					else if (i > j)
						result[i][j] = matrix[i][j];
				}
			}
			return result;

		}

		static double[][] ExtractUpper(double[][] matrix)
		{
			// upper part of a Doolittle decomposition (0.0s in the strictly lower part)
			int rows = matrix.Length; int cols = matrix[0].Length;
			double[][] result = MatrixCreate(rows, cols);
			for (int i = 0; i < rows; ++i)
			{
				for (int j = 0; j < cols; ++j)
				{
					if (i <= j)
						result[i][j] = matrix[i][j];
				}
			}
			return result;
		}

		static double[][] PermArrayToMatrix(int[] perm)
		{
			// convert Doolittle perm array to corresponding perm matrix
			int n = perm.Length;
			double[][] result = MatrixCreate(n, n);
			for (int i = 0; i < n; ++i)
				result[i][perm[i]] = 1.0;
			return result;
		}

		static double[][] UnPermute(double[][] luProduct, int[] perm)
		{
			// unpermute product of Doolittle lower * upper matrix according to perm[]
			// no real use except to demo LU decomposition, or for consistency testing
			double[][] result = MatrixDuplicate(luProduct);

			int[] unperm = new int[perm.Length];
			for (int i = 0; i < perm.Length; ++i)
				unperm[perm[i]] = i;

			for (int r = 0; r < luProduct.Length; ++r)
				result[r] = luProduct[unperm[r]];

			return result;
		}

		static string VectorAsString(double[] vector)
		{
			string s = "";
			for (int i = 0; i < vector.Length; ++i)
				s += vector[i].ToString("F3").PadLeft(8) + Environment.NewLine;
			s += Environment.NewLine;
			return s;
		}

		static string VectorAsString(int[] vector)
		{
			string s = "";
			for (int i = 0; i < vector.Length; ++i)
				s += vector[i].ToString().PadLeft(2) + " ";
			s += Environment.NewLine;
			return s;
		}
	}
}
