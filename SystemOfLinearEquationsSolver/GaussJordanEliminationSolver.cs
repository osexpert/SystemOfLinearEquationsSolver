using MathNet.Numerics;
using System;
using System.Security.Cryptography;
using System.Threading;
using System.Timers;

namespace SystemOfLinearEquationsSolver
{
	

	public class GaussJordanEliminationSolver : ISolver
	{

		public static readonly GaussJordanEliminationSolver Instance = new GaussJordanEliminationSolver();

        private GaussJordanEliminationSolver()
        {
            
        }


		/// <summary>
		/// 
		/// Example with 2 unknowns:
		/// 1x + 2y = 4
		/// 2x + 3y = 5
		/// Coefficients will be:
		/// new double[,] 
		/// {
		///		{1, 2, 4},
		///		{2, 3, 5}
		///	};
		///	The result will be: { -2, 3 }
		///	x = -2, y = 3
		///	
		/// Gaussian elimination with partial pivoting to transform the input matrix into reduced row-echelon form. 
		/// The solutions of the system of linear equations are then extracted from the last column of the reduced matrix.
		/// Keep in mind that while Gaussian elimination is a powerful method for solving systems of linear equations, 
		/// it might not handle all possible cases.In some scenarios, the matrix might be ill-conditioned or singular,
		/// leading to numerical instability or division by zero issues.Additionally, this method might not be the most 
		/// efficient for large matrices, as it has a cubic time complexity.
		/// 
		/// Changed back from jagged to non-jagged because swap rows was so expensive.
		/// </summary>
		/// <param name="_coefficients"></param>
		/// <returns></returns>
		public double[] SolveEquations(double[,] coefficients, int? resultRoundDigits = null)
		{
			return SolveEquations(SolverCommon.GetMultiFromJaggedCopy(coefficients), resultRoundDigits);
		}

		double[] SolveEquations(double[][] coefficientsCopy, int? resultRoundDigits = null)
		{ 
			int numRows = coefficientsCopy.Length;
			int numCols = coefficientsCopy[0].Length;
			for (int pivot = 0; pivot < numRows; pivot++)
			{
				// Find the row with the largest absolute pivot value below the current pivot row
				int maxRow = pivot;
				double maxPivotValue = Math.Abs(coefficientsCopy[pivot][pivot]);
				for (int row = pivot + 1; row < numRows; row++)
				{
					double currentPivotValue = Math.Abs(coefficientsCopy[row][pivot]);
					if (currentPivotValue > maxPivotValue)
					{
						maxPivotValue = currentPivotValue;
						maxRow = row;
					}
				}

				// Swap the current pivot row with the row having the largest pivot value
				if (maxRow != pivot)
				{
					// Swap row
					var temp = coefficientsCopy[pivot];
					coefficientsCopy[pivot] = coefficientsCopy[maxRow];
					coefficientsCopy[maxRow] = temp;
				}

				// Normalize the pivot row by dividing all elements by the pivot value
				double pivotValue = coefficientsCopy[pivot][pivot];
				if (pivotValue == 0d)
					throw new SolverException("pivotValue is 0, would divide by zero");

				for (int col = 0; col < numCols; col++)
				{
					coefficientsCopy[pivot][col] /= pivotValue;
				}

				// Use the pivot row to eliminate the coefficients in other rows
				for (int row = 0; row < numRows; row++)
				{
					if (row != pivot)
					{
						double factor = coefficientsCopy[row][pivot];
						for (int col = 0; col < numCols; col++)
						{
							// Subtract the pivot row multiplied by the factor from the current row
							coefficientsCopy[row][col] -= factor * coefficientsCopy[pivot][col];
						}
					}
				}
			}

			// Extract the solutions from the last column of the reduced matrix
			double[] solution = new double[numRows];
			for (int i = 0; i < numRows; i++)
			{
				if (resultRoundDigits != null)
					solution[i] = Math.Round(coefficientsCopy[i][numRows], resultRoundDigits.Value);
				else
					solution[i] = coefficientsCopy[i][numRows];
			}

			return solution;
		}


	}
}
