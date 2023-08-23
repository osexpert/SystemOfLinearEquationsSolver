using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SystemOfLinearEquationsSolver
{
	public class SubstitutionSolver
	{
		/// <summary>
		/// Solve fast path where one equation has only one unknown and can be solved directly, and all other equations can be solved with direct substitution.
		/// If not, do some other slow path solver.
		/// 
		/// coeff: Ax=B
		/// </summary>
		/// <param name="_coeff"></param>
		/// <returns></returns>
		public static double[] SolveHappyPathOrFailFast(double[,] coeff)
		{
			var rows = coeff.GetLength(0);
			var cols = coeff.GetLength(1) - 1; // 0=x, 1=y, etc,

			int numColsSolved = 0;
			var rowSolved = new bool[cols];
			var solutions = new double[cols];

			bool solved_anything_in_this_pass = false;
			do
			{
				solved_anything_in_this_pass = false;

				for (int i = 0; i < rows; i++)
				{
					if (rowSolved[i])
						continue;

					// find cols without answers. these are cols with coeff != 0 and cols that not already have answers.
					double sum = 0;

					int? single_col_needs_answer = null;

					for (int j = 0; j < cols; j++)
					{
						
						if (coeff[i, j] == 0d)
						{
							// no unknown in this col, ignore
						}
						else if (rowSolved[j])
						{
							sum += coeff[i, j] * solutions[j];
						}
						else
						{
							// col needs an answer
							if (single_col_needs_answer == null)
							{
								single_col_needs_answer = j;
							}
							else
							{
								single_col_needs_answer = null;
								break; // more than one col needs answer, skip
							}
						}
					}

					// have a single col that need asweer
					if (single_col_needs_answer != null)
					{
						var contantEquals = coeff[i, cols];
						var solution_to_col = (contantEquals - sum) / coeff[i, single_col_needs_answer.Value];

						if (rowSolved[single_col_needs_answer.Value])
							throw new Exception("Bug, already solved");

						solutions[single_col_needs_answer.Value] = solution_to_col;

						solved_anything_in_this_pass = true;
						rowSolved[i] = true;
						numColsSolved++;

						if (numColsSolved == cols)
							break;
					}

				}

			} while (solved_anything_in_this_pass && numColsSolved < cols);

			// if we have all solutions we are golden, else return null
			if (numColsSolved == cols)
				return solutions;
			else
				return null;
		}
	}
}
