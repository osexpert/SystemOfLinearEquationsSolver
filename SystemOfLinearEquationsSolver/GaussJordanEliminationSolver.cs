using System;
using System.Threading;

namespace SystemOfLinearEquationsSolver
{
	

	public class GaussJordanEliminationSolver : ISolver
	{

		public static readonly GaussJordanEliminationSolver Instance = new GaussJordanEliminationSolver();

        private GaussJordanEliminationSolver()
        {
            
        }

        /// <summary>
        /// Gauss-Jordan elimination with back-substitution?
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
		/// Works fine, but a bit slow with many unknows, it does not scale well.
        /// </summary>
        /// <param name="_coefficients"></param>
        /// <param name="numEquations"></param>
        /// <returns></returns>
        public double[] SolveEquations(double[,] coefficients, int? resultRoundDigits = null)
		{
			var coefficientsCopy = Copy(coefficients);
			int numEquations = coefficientsCopy.GetLength(0);

			for (int pivot = 0; pivot < numEquations; pivot++)
			{
				double pivotValue = coefficientsCopy[pivot, pivot];


				// Normalize the pivot row
				for (int col = 0; col <= numEquations; col++)
				{
					if (pivotValue == 0d)
						throw new SolverException("pivotValue is 0, would divide by zero");

					coefficientsCopy[pivot, col] /= pivotValue;
				}


				// Eliminate other rows
				for (int row = 0; row < numEquations; row++)
				{
					if (row != pivot)
					{
						double factor = coefficientsCopy[row, pivot];
						for (int col = 0; col <= numEquations; col++)
						{
							coefficientsCopy[row, col] -= factor * coefficientsCopy[pivot, col];
						}
					}
				}
			}

			// Extract solutions (back-substitution)
			double[] solution = new double[numEquations];
			for (int i = 0; i < numEquations; i++)
			{
				if (resultRoundDigits != null)
					solution[i] = Math.Round(coefficientsCopy[i, numEquations], resultRoundDigits.Value);
				else
					solution[i] = coefficientsCopy[i, numEquations];
			}

			return solution;
		}

		private double[,] Copy(double[,] coefficients)
		{
			// TODO: can use Clone()?
			var a = coefficients.GetLength(0);
			var b = coefficients.GetLength(1);

			var res = new double[a, b];

			for (int i = 0; i < a; i++)
			{
				for (int j = 0; j < b; j++)
				{
					res[i, j] = coefficients[i, j];
				}
			}

			return res;
		}

	}
}
