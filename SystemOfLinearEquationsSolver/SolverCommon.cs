using System;
using System.Collections.Generic;
using System.Text;

namespace SystemOfLinearEquationsSolver
{
	public class SolverCommon
	{

		public static (double[,], double[]) GetTestData(int unknowns, bool firstRowOnlyOneUnknown)
		{
			Random r = new Random();

			// generate unknowns (x, y, z)
			double[] unk = new double[unknowns];
			for (int i = 0; i < unknowns; i++)
			{
				var d = r.NextDouble();
				unk[i] = d > 0.5 ? d * 100 : -d * 100;
			}

			double[,] res = new double[unknowns, unknowns + 1];
			//double[,] res2 = new double[unknowns, unknowns];
			for (int i = 0; i < unknowns; i++)
			{
				double sum = 0d;
				for (int j = 0; j < unknowns; j++)
				{
					// diagonal ans straight down først column
					if (i == j || j == 0 || i == 0 && j == 1 && !firstRowOnlyOneUnknown)
					{
						var d = r.NextDouble();
						//var d2 = r.NextDouble();
						res[i, j] = d > 0.5 ? d * 100 : -d * 100; // generate value of unknown (Nx)
						//res2[i, j] = d2 > 0.5 ? d2 * 100 : -d2 * 100;
						// multiply these two and get the anwwer
						sum += res[i, j] * unk[j];
					}

					
				}

				res[i, unknowns] = sum;
			}

			return (res, unk);
		}

		public static bool CloseToZero(double d, int doublesToCheck)
		{

			return Math.Abs(d) < Math.Pow(10, -doublesToCheck);
		}

		/// <summary>
		/// Insert result in the equations and sum it up and compare. Any deviation is returned.
		/// Ideally, all deviations should be 0, but sometimes there can be fractional errors.
		/// </summary>
		/// <param name="coefficients"></param>
		/// <param name="result"></param>
		/// <returns></returns>
		public static double[] GetDeviations(double[,] coefficients, double[] result)
		{
			double[] deviations = new double[result.Length];

			var dim0Length = coefficients.GetLength(0);
			var dim1Length = coefficients.GetLength(1);

			for (int i = 0; i < dim0Length; i++)
			{
				double sum = 0d;

				for (int j = 0; j < dim1Length; j++)
				{
					if (j == dim1Length - 1)
					{
						sum -= coefficients[i, j];
					}
					else
					{
						sum += coefficients[i, j] * result[j];
					}

				}

				deviations[i] = sum;
			}

			return deviations;
		}

	}

	public interface ISolver
	{
		/// <summary>
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
		/// Often result will contain small fractional errors, because of division of double is a losse calculation.
		/// Rounding the result to 14-15 digits will often suffice to remove the fractional errors.
		/// </summary>
		/// <param name="coefficients"></param>
		/// <param name="resultRoundDigits"></param>
		/// <returns></returns>
		double[] SolveEquations(double[,] coefficients, int? resultRoundDigits = null);
	}

	public class SolverException : ApplicationException
	{
		public SolverException(string msg) : base(msg)
		{

		}
	}
}
