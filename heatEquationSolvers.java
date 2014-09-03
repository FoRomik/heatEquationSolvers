import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Scanner;

public class HW3_Driver
{
	private static final double MIN = 1e-8;  // Minimum error

	// Initial Gaussin condition function with or without noise
	public static double initialCondition(double x, double y, double z, double x_mean, double y_mean, double z_mean, boolean noiseBoolean)
	{
		Random rand = new Random();
		double initial = Math.exp(-Math.pow(x - x_mean, 2) - Math.pow(y - y_mean , 2) - Math.pow(z - z_mean, 2));

		if (noiseBoolean)
			initial +=  0.00001 * (double)(rand.nextInt(100) - 50);

		return initial;
	}

	// Source term function
	public static double source(double i, double j, double k)
	{	
		return 0.001 * (i + j + k);
	}

	// Convert 3D matrix to vector
	public static double[] threeDtoVector(double[][][] T, int nx, int ny, int nz)
	{
		double[] TVector = new double[(nx - 1) * (ny - 1) * (nz - 1)];

		int vectorIndex = 0;

		for(int i = 1; i < nx; i ++)
		{
			for (int j = 1; j < ny; j++)
			{
				for (int k = 1; k < nz; k++)
				{
					TVector[vectorIndex] = T[i][j][k];
					vectorIndex++;
				}
			}
		}
		return TVector;
	}


	// Given indexes in 3D form, return a 1D vector index
	public static int index3DtoVector(int i, int j, int k, int nx, int ny, int nz)
	{
		int ret = (i - 1) * (ny - 1) * (nz - 1) + (j - 1) * (nz - 1) + k - 1;
		return ret;
	}

	// Convert a 1D T vector to 3D form
	public static double[][][] vectorTo3D(double[] TVector, int nx, int ny, int nz, double boundConst, boolean dirichlet)
	{
		double[][][] T = new double[nx + 1][ny + 1][nz + 1];

		int vectorIndex = 0;

		//Interior points
		for(int i = 1; i < nx; i++)
		{
			for(int j = 1; j < ny; j++)
			{
				for(int k = 1; k < nz; k++)
				{
					T[i][j][k] = TVector[vectorIndex];
					vectorIndex++;
				}
			}
		}


		// Boundary points
		for(int i = 0; i <= nx; i++)
		{
			for(int j = 0; j <= ny; j++)
			{
				for(int k = 0; k <= nz; k++)
				{
					if(i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz)
					{	
						if(dirichlet)  // Constant boundary value
						{
							T[i][j][k] = boundConst;
						}
						else  // Periodic boundary value
						{
							int p, q, r;
							p = i;
							q = j;
							r = k;

							if(p == 0)
								p = nx - 1;
							if(p == nx)
								p = 1;
							if(q == 0)
								q = ny - 1;
							if(q == ny)
								q = 1;
							if(r ==0)
								r = nz - 1;
							if(r == nz)
								r = 1;

							T[i][j][k] = T[p][q][r];
						}
					}
				}
			}
		}
		return T;
	}

	// Initialize a Temperature matrix in time 0
	public static double[][][] initialization(double dx, int nx, int ny, int nz, boolean dirichlet, double boundConst, boolean noiseBoolean)
	{
		double x_mean = dx * nx / 2;
		double y_mean = dx * ny / 2;
		double z_mean = dx * nz / 2;

		double[][][] T_init = new double[nx + 1][ny + 1][nz + 1];


		// Initializing interior condition
		for(int i = 1; i < nx; i++)
		{
			for(int j = 1; j < ny; j++)
			{
				for(int k = 1; k < nz; k++)
				{
					T_init[i][j][k] = initialCondition(i * dx, j * dx, k * dx, x_mean, y_mean, z_mean, noiseBoolean);

				}
			}
		}

		// Initializing boundary condition
		for(int i = 0; i <= nx; i++)
		{
			for(int j = 0; j <= ny; j++)
			{
				for(int k = 0; k <= nz; k++)
				{
					if(i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz)
					{	
						if(dirichlet)
						{
							T_init[i][j][k] = boundConst;
						}
						else
						{
							int p, q, r;
							p = i;
							q = j;
							r = k;

							if(p == 0)
								p = nx - 1;
							if(p == nx)
								p = 1;
							if(q == 0)
								q = ny - 1;
							if(q == ny)
								q = 1;
							if(r ==0)
								r = nz - 1;
							if(r == nz)
								r = 1;

							T_init[i][j][k] = T_init[p][q][r];
						}
					}
				}
			}

		}

		return T_init;
	}



	public static void CrankNicolsonMethod(double[][][] T, double dt, double dx, double alpha, 
			int nx, int ny, int nz, int nt, boolean dirichlet, double boundConst, double lambda, boolean sourceBoolean, int method)
	{
		double[] TVector = new double[(nx - 1) * (ny - 1) * (nz - 1)];
		double[] b = new double[(nx - 1) * (ny - 1) * (nz - 1)];
		//double[][] A = new double[(nx - 1) * (ny - 1) * (nz - 1)][(nx - 1) * (ny - 1) * (nz - 1)];
		int[][]A = new int[(nx - 1) * (ny - 1) * (nz - 1)][6];   // A[i][j] is the column index of non-zero elements on row i of the original A matrix  


		switch (method)
		{
		case 1:
			System.out.println("\nJacobi Method:");
			break;
		case 2:
			System.out.println("\nGauss-Seidel Method:");
			break;
		case 3:
			System.out.println("\nSOR Method:");
			break;
		default:
			break;
		}

		TVector = threeDtoVector(T, nx, ny, nz);

		int numBoundPoints = 0;

		for(int rowIndex = 0; rowIndex < (nx - 1) * (ny - 1) * (nz - 1); rowIndex++)
		{
			int p = rowIndex / ((ny - 1) * (nz - 1)) + 1;
			int q = (rowIndex % ((nz - 1) * (ny - 1))) / (nz - 1) + 1;
			int r = (rowIndex % ((nz - 1) * (ny - 1))) % (nz - 1) + 1;

			if (p > 1)
				A[rowIndex][0] = index3DtoVector(p - 1, q, r, nx, ny, nz);
			else
				A[rowIndex][0] = 0;

			if (p < nx - 1)
				A[rowIndex][1] = index3DtoVector(p + 1, q, r, nx, ny, nz);
			else
				A[rowIndex][1] = 0;

			if (q > 1)
				A[rowIndex][2] = index3DtoVector(p, q - 1, r, nx, ny, nz);
			else
				A[rowIndex][2] = 0; 

			if (q < ny - 1)
				A[rowIndex][3] = index3DtoVector(p, q + 1, r, nx, ny, nz);
			else 
				A[rowIndex][3] = 0; 

			if (r > 1)
				A[rowIndex][4] = index3DtoVector(p, q, r - 1, nx, ny, nz);
			else
				A[rowIndex][4] = 0; 

			if (r < nz - 1)
				A[rowIndex][5] = index3DtoVector(p, q, r + 1, nx, ny, nz);
			else
				A[rowIndex][5] = 0; 			
		}



		int theta = 0;
		System.out.format("%5d", theta);
		for(int i = 0; i <= nx; i++)
		{
			System.out.format("%10.6f ", T[i][(int)(ny / 2)][(int)(nz / 2)]);
		}
		System.out.println();

		// Updating each time step
		for (theta = 1; theta <= nt; theta++)
		{
			System.out.format("%5d", theta);

			int rowIndex = 0;

			for(int i = 1; i < nx; i++)
			{
				for(int j = 1; j < ny; j++)
				{
					for(int k = 1; k < nz; k++)
					{
						b[rowIndex] = 2 * T[i][j][k] + lambda * ( T[i - 1][j][k] + T[i + 1][j][k] + T[i][j - 1][k] 
								+ T[i][j + 1][k] + T[i][j][k - 1] + T[i][j][k + 1] - 6 * T[i][j][k]); 

						if (dirichlet)
						{
							numBoundPoints = 0;
							if (i == 1 || i == nx - 1)
								numBoundPoints++;
							if (j == 1 || j == ny - 1)
								numBoundPoints++;
							if (k == 1 || k == nz - 1)
								numBoundPoints++;
							b[rowIndex] = b[rowIndex]+ numBoundPoints * boundConst * lambda;
						}

						if (sourceBoolean)
							T[i][j][k] += source(i,j,k);

						rowIndex++;
					}
				}
			}

			//TVector = GaussJordanSolve.solve(A, b);
			TVector = iterativeSolve(A, b, method, lambda);

			T = vectorTo3D(TVector, nx, ny, nz, boundConst, dirichlet);

			// Show
			for(int i = 0; i <= nx; i++)
			{
				System.out.format("%10.6f ", T[i][(int)(ny / 2)][(int)(nz / 2)]);
			}
			System.out.println();
		}
	}


	public static double[] iterativeSolve(int[][] A, double[] b, int method, double lambda)
	{
		int iterMax = 1000;

		int length = A.length;

		double[] T_old = new double[length];

		double[] T_new = new double[length];

		for (int i = 0; i < iterMax; i++)
		{
			switch (method)
			{
			case 1:
				T_new = jacobi (A, b, T_old, lambda);
				break;
			case 2:
				T_new = gaussSeidel (A, b, T_old, lambda);
				break;
			case 3:
				T_new = sor (A, b, T_old, lambda);
				break;
			default:
				break;
			}


			double abs_sum = 0;
			for (int j = 0; j < length; j++)
			{
				abs_sum = Math.abs(T_new[j] - T_old[j]);
			}
			if (abs_sum / length < MIN)
				break;


			for (int j = 0; j < length; j++)
				T_old[j] = T_new[j];			
		}

		return T_new;
	}


	public static double[] jacobi (int[][] A, double[] b, double[] T_old, double lambda)
	{
		int length = T_old.length;

		double[] T_new = new double[length];

		for(int i = 0; i < length; i++)
		{
			double sum = 0;
			for(int j = 0; j < 6; j++)
			{
				if (A[i][j] > 0)
					sum += T_old[A[i][j]] * (- lambda); 

			}

			T_new[i] = (b[i] - sum) / (2 + 6 * lambda);

		}
		return T_new;
	}


	public static double[] gaussSeidel (int[][] A, double[] b, double[] T_old, double lambda)
	{
		int length = T_old.length;

		double[] T_new = new double[length];

		for(int i = 0; i < length; i++)
		{
			double sum1 = 0;
			double sum2 = 0;


			for(int j = 0; j < 6; j++)
			{
				if (A[i][j] > 0)
				{
					if(A[i][j] < i)
						sum1 += T_new[A[i][j]] * (- lambda); 
					else if (A[i][j] > i)
						sum2 += T_old[A[i][j]] * (- lambda);
				}

			}

			T_new[i] = (b[i] - sum1 - sum2) / (2 + 6 * lambda);

		}
		return T_new;
	}


	public static double[] sor (int[][] A, double[] b, double[] T_old, double lambda)
	{
		double omega = 1.75;
		int length = T_old.length;

		double[] T_new = new double[length];

		for(int i = 0; i < length; i++)
		{
			double sum1 = 0;
			double sum2 = 0;


			for(int j = 0; j < 6; j++)
			{
				if (A[i][j] > 0)
				{
					if(A[i][j] < i)
						sum1 += T_new[A[i][j]] * (- lambda); 
					else if (A[i][j] > i)
						sum2 += T_old[A[i][j]] * (- lambda);
				}

			}

			T_new[i] =  (1 - omega) * T_old[i] + omega * (b[i] - sum1 - sum2) / (2 + 6 * lambda);

		}
		return T_new;
	}

	public static void main(String[] args) 
	{
		double dt, dx, alpha, lambda, boundConst;
		int nx, ny, nz, nt;
		boolean dirichlet, sourceBoolean, noiseBoolean;
		String userInput;

		double[][][] T_init1;
		double[][][] T_init2;
		double[][][] T_init3;

		dirichlet = true;
		sourceBoolean = false;
		noiseBoolean = false;

		dx = 0.1;
		dt = 0.0025;
		nx = 10;
		ny = 10;
		nz = 10;
		nt = 10;
		alpha = 0.5;
		boundConst = 0;

		Scanner input = new Scanner(System.in);

		do{
			System.out.println("Runtime performace evaluation for all methods? (Y/N)");
			userInput = input.next(); 
		} while (!(userInput.toUpperCase().equals("Y") || userInput.toUpperCase().equals("N")));
		if (userInput.toUpperCase().equals("Y"))
		{
			try 
			{
				File file = new File("runtimeEvaluation.txt");

				if (!file.exists()) 
				{
					file.createNewFile();
				}

				FileWriter fw = new FileWriter(file.getAbsoluteFile());
				BufferedWriter bw = new BufferedWriter(fw);

				long[] duration = new long[65];
				long start, end;
				start = System.currentTimeMillis();

				nt = 1;
				nx = 1;

				for (int i = 0; i < duration.length; i++)
				{
					nx += 2;
					ny = nz = nx;
					lambda = alpha * dt / (Math.pow(dx, 2));

					bw.write(nx + " ");

					T_init1= initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);
					T_init2= initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);
					T_init3= initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);

					start = System.currentTimeMillis();
					CrankNicolsonMethod(T_init1, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 1);
					end = System.currentTimeMillis();
					duration[i] = end - start;
					System.out.println("nx =  " + nx + "; duration: " + duration[i] + "ms");
					bw.write(duration[i] + " ");

					start = System.currentTimeMillis();
					CrankNicolsonMethod(T_init2, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 2);
					end = System.currentTimeMillis();
					duration[i] = end - start;
					System.out.println("nx =  " + nx + "; duration: " + duration[i] + "ms");
					bw.write(duration[i] + " ");

					start = System.currentTimeMillis();
					CrankNicolsonMethod(T_init3, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 3);
					end = System.currentTimeMillis();
					duration[i] = end - start;
					System.out.println("nx =  " + nx + "; duration: " + duration[i] + "ms");
					bw.write(duration[i] + "\n");
				}

				bw.close();
				System.out.println("Done.");
			}catch (IOException e) {
				e.printStackTrace();
			}
		}
		else
		{
			System.out.println("dx = " + dx + "\ndt = " + dt + "\nnx = " + nx + "\nny = " + ny + "\nnz " + nz + "\nnt = " + nt + 
					"\nalpha = " + alpha + "\nboundConst " + boundConst + "\ndirichlet = " + dirichlet + "\nsourceBoolean = " 
					+ sourceBoolean + "\nnoiseBoolean = " + noiseBoolean);


			lambda = alpha * dt / (Math.pow(dx, 2));

			System.out.format("lambda = %5.6f\n", lambda);

			T_init1 = new double[nx + 1][ny + 1][nz + 1];
			T_init2 = new double[nx + 1][ny + 1][nz + 1];
			T_init3 = new double[nx + 1][ny + 1][nz + 1];

			T_init1 = initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);
			T_init2 = initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);
			T_init3 = initialization(dx, nx, ny, nz, dirichlet, boundConst, noiseBoolean);

			CrankNicolsonMethod(T_init1, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 1);
			CrankNicolsonMethod(T_init2, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 2);
			CrankNicolsonMethod(T_init3, dt, dx, alpha, nx, ny, nz, nt, dirichlet, boundConst, lambda, sourceBoolean, 3);
		}
	}
}