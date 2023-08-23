

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Diagnostics;
using SystemOfLinearEquationsSolver;

namespace test_net6
{
	internal class Program
	{
		static void Main(string[] args)
		{
			//subsolver.Maine(null);

			// if firstRowOnlyOneUnknown: true, we have a happy path
			var t = SolverCommon.GetTestData(1000, firstRowOnlyOneUnknown: false);

			var s = Stopwatch.StartNew();

			//10 sec on 1000 unk. same with 3unk.
			//var res = GaussJordanEliminationSolver.Instance.SolveEquations(t.Item1);

			//2.5 sec on 1000 unk.
			//3sec with 2unk per row
			//var res = MatrixSolver.Instance.SolveEquations(t.Item1);

			// 110ms for 1000 unk... KE??? Yes, when row1 only has 1 unk it can solve it fast...happy path optimization.
			// 22sek, with 2 unk per row. It could have failed over to matrix for complex solutons (it does, but it is slow and recursive)
			//						var sol = new ClassicSolver(GetMultiFromJagged(t.Item1));
			//				var res = sol.SolveEquation();

			// approx 65ms
			var res = SubstitutionSolver.SolveHappyPathOrFailFast(t.Item1);
			if (res == null)
			{
				//				var ab = SplitIntoAx_B(t.Item1);
				//			res = ab.Item1.Solve(ab.Item2).ToArray(); // if (ColumnCount == RowCount) LU else QR...

				// Seem to use same time as MathNet, but I guess it make sense.
				res = MatrixSolver_LU.Instance.SolveEquations(t.Item1);
			}


			//3,2sek with firstRowOnlyOneUnk: false, same with true

			//3,2sek with firstRowOnlyOneUnk: false, same with true
			//var res = ab.Item1.LU().Solve(ab.Item2).ToArray();

			//2,4 - 2,6. The fastest?
			//var res = ab.Item1.QR().Solve(ab.Item2).ToArray();

			//15sek, firstRowOnlyOneUnk: false
			//var res = ab.Item1.Svd().Solve(ab.Item2).ToArray();

			//2.2sek, firstRowOnlyOneUnk: false
			//var res = ab.Item1.GramSchmidt().Solve(ab.Item2).ToArray();

			//did not work
			//var res = ab.Item1.Evd().Solve(ab.Item2).ToArray();

			s.Stop();

			Console.WriteLine("Time: ms " + s.ElapsedMilliseconds);

			var devs = SolverCommon.GetDeviations(t.Item1, res);
			foreach (var d in devs)
				if (!SolverCommon.CloseToZero(d, 10))
				{
				}
		}

		private static (MathNet.Numerics.LinearAlgebra.Double.DenseMatrix, Vector<double>) SplitIntoAx_B(double[,] item1)
		{
			double[] conta = new double[item1.GetLength(0)];
			double[,] ama = new double[item1.GetLength(0), item1.GetLength(1) - 1];

			for (int i = 0; i < item1.GetLength(0); i++)
			{
				for (int j = 0; j < item1.GetLength(1) - 1; j++)
				{
					ama[i, j] = item1[i, j];
				}
				conta[i] = item1[i, item1.GetLength(1) - 1];
			}

			var m1 = DenseMatrix.OfArray(ama);

			//double[] loadVector = { 0, 0, 0, 5000, 0, 0, 0, 0, 0 };
			var vectorB = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(conta);

			return (m1, vectorB);
		}

		private static double[][] GetMultiFromJagged(double[,] arr)
		{
			var dim1len = arr.GetLength(1);
			double[][] res = new double[arr.GetLength(0)][];
			for (int i = 0; i < res.Length; i++)
			{
				res[i] = new double[dim1len];
				for (int j = 0; j < dim1len; j++)
				{
					res[i][j] = arr[i, j];
					if (j == dim1len - 1)
						res[i][j] = -res[i][j];
				}
			}
			return res;
		}
	}
}