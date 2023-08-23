

using SystemOfLinearEquationsSolver;

namespace UnitTests
{
	class Test : IDisposable
	{
		private bool disposedValue;

		protected virtual void Dispose(bool disposing)
		{
			if (!disposedValue)
			{
				if (disposing)
				{
					// TODO: dispose managed state (managed objects)
				}

				// TODO: free unmanaged resources (unmanaged objects) and override finalizer
				// TODO: set large fields to null
				disposedValue = true;
			}
		}

		// // TODO: override finalizer only if 'Dispose(bool disposing)' has code to free unmanaged resources
		// ~Test()
		// {
		//     // Do not change this code. Put cleanup code in 'Dispose(bool disposing)' method
		//     Dispose(disposing: false);
		// }

		public void Dispose()
		{
			// Do not change this code. Put cleanup code in 'Dispose(bool disposing)' method
			Dispose(disposing: true);
			GC.SuppressFinalize(this);
		}
	}

	public class UnitTest1
	{
		[Theory]
		[InlineData(enMethod.GaussJordan, 14)]
		[InlineData(enMethod.Matrix_LU, 14)]
		public void Test1(enMethod m, int? roundDigits = null)
		{
			double[,] coefficients = {
					{ 1, 1, 0, 7 },
					{ -2, -1, 2, 0 },
					{ 3, -2, -2, -9 }
				};

			var res = GetFromMeth(m).SolveEquations(coefficients, roundDigits);
			Assert.Equal(new double[] { 3, 4, 5 }, res);
			Verify(coefficients, res);
		}

		private void Verify(double[,] coefficients, double[] results)
		{
			var res = SolverCommon.GetDeviations(coefficients, results);
			foreach (var d in res)
			{
				Assert.Equal(0d, d);
			}
		}

		/// <summary>
		/// In this test, the Matrix solver seems to be much slower than the GaussJordan solver...
		/// </summary>
		/// <param name="m"></param>
		/// <param name="roundDigits"></param>
		[Theory]
		[InlineData(enMethod.GaussJordan, 14)]
		[InlineData(enMethod.Matrix_LU, 14)]
		public void Test2(enMethod m, int? roundDigits = null)
		{

			double[,] coefficients = {
				{ 3, 2, -1, 1 },
				{ -1, 0.5, -1, 0 },
				{ 2, -2, 4, -2 }
			};

			var res = GetFromMeth(m).SolveEquations(coefficients, roundDigits);
			Assert.Equal(new double[] { 1, -2, -2 }, res);
			//1, -2, -2 (needed rounding to get it right)
			Verify(coefficients, res);
		}

		[Theory]
		[InlineData(enMethod.GaussJordan, 15)]
		[InlineData(enMethod.Matrix_LU, 15)]
		public void Test3(enMethod m, int? roundDigits = null)
		{
			double[,] coefficients = {
				{ 2, 1, -1, 8 },
				{ -3, -1, 2, -11 },
				{ -2, 1, 2, -3 }
			};

			var res = GetFromMeth(m).SolveEquations(coefficients, roundDigits);
			Assert.Equal(new double[] { 2, 3, -1 }, res);
			Verify(coefficients, res);
		}

		[Theory]
		[InlineData(enMethod.GaussJordan)]
		[InlineData(enMethod.Matrix_LU)]
		public void Test4(enMethod m)
		{
			double[,] coefficients = {
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, 0 }
			};

			Assert.Throws<SolverException>(() =>
			{
				var res = GetFromMeth(m).SolveEquations(coefficients);
			});


		}

		[Theory]
		[InlineData(enMethod.GaussJordan)]
		[InlineData(enMethod.Matrix_LU)]
		public void Test5(enMethod m)
		{
			double[,] coefficients = {
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, double.NaN }
			};

			Assert.Throws<SolverException>(() =>
			{
				var res = GetFromMeth(m).SolveEquations(coefficients);
			});


		}

		[Theory]
		[InlineData(enMethod.GaussJordan)]
		[InlineData(enMethod.Matrix_LU)]
		public void Test6(enMethod m)
		{
			double[,] coefficients = {
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, double.PositiveInfinity }
			};

			Assert.Throws<SolverException>(() =>
			{
				var res = GetFromMeth(m).SolveEquations(coefficients);
			});


		}

		[Theory]
		[InlineData(enMethod.GaussJordan)]
		[InlineData(enMethod.Matrix_LU)]
		public void Test7(enMethod m)
		{
			double[,] coefficients = {
				{ 1, 1, 1, 1 },
				{ 1, 1, 1, double.NegativeInfinity }
			};

			Assert.Throws<SolverException>(() =>
			{
				var res = GetFromMeth(m).SolveEquations(coefficients);
			});


		}

		[Theory]
		[InlineData(enMethod.GaussJordan)]
		[InlineData(enMethod.Matrix_LU)]
		public void Test8(enMethod m)
		{
			double[,] coefficients = {
				{1, 2, 4},
				{2, 3, 5}
			};

			var res = GetFromMeth(m).SolveEquations(coefficients);
			Assert.Equal(new double[] { -2, 3 }, res);
			Verify(coefficients, res);
		}



		public static ISolver GetFromMeth(enMethod m)
		{
			switch (m)
			{
				case enMethod.GaussJordan: return GaussJordanEliminationSolver.Instance;
				case enMethod.Matrix_LU: return MatrixSolver.Instance;
				default:
					throw new NotImplementedException();
			}
		}

		public enum enMethod
		{
			GaussJordan = 1,
			Matrix_LU = 2,
		}
	}
}