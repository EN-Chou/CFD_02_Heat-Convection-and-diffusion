// CDS currently fails at pe>1
#include "header.h"

//Parameters
double m_grid[N][N] = { 0 };
double h = 1 / ((double)N - 1);
double u = 1, v = 1, rho = 1;
double gamma = (double)u / (Pe * (N - 1));
double m_rate = rho * u * h;
double D_e = gamma, D_w = gamma, D_n = gamma, D_s = gamma;
double r_e = 0, r_w = 0, r_n = 0, r_s = 0;
double phi_e = 0, phi_w = 0, phi_n = 0, phi_s = 0;
double A_e, A_w, A_n, A_s, A_ww, A_ss, A_p;
double time2;
int i = 0, j = 0;
int iterations = 0;
double residual = 1;
double residual_after = 1;
double temp = 0;
int partition[10];


int main(int argc, char* argv[]) {
	int rank;
	int size;
	MPI_Status istat[8];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*Implementation*/
	//CPU distribution
	for (i = 0; i <= size; i++) {
		partition[i] = (N - 1) / size * i;
	}

	//Initialization
	for (i = 0; i < N; i++) {
		m_grid[i][0] = 1.0;
		m_grid[i][N - 1] = 0.0;
		m_grid[N - 1][i] = 0.0;
		m_grid[0][i] = 1.0;
	}

	//main
	double time1 = MPI_Wtime();
	double temp2;
	while (residual > tol) {
		residual = 0;
		iterations++;

		/*2nd outer layer by upwind*/
		phi_generate(0);
		for (i = partition[rank]; i < partition[rank + 1]; i++) {
			if (i == N - 2) {
				for (j = 1; j < N - 1; j++) {
					temp = m_grid[i][j];
					m_grid[i][j] = (m_grid[i][j + 1] * A_e + m_grid[i][j - 1] * A_w + m_grid[i - 1][j] * A_n + m_grid[i + 1][j] * A_s) / A_p;
					residual = residual + abs(m_grid[i][j] - temp);
				}
			}

			if (i != 0) {
				j = 1;
				temp = m_grid[i][j];
				m_grid[i][j] = (m_grid[i][j + 1] * A_e + m_grid[i][j - 1] * A_w + m_grid[i - 1][j] * A_n + m_grid[i + 1][j] * A_s) / A_p;
				residual = residual + abs(m_grid[i][j] - temp);
			}
		}

		//inner
		for (i = partition[rank]; i < partition[rank + 1]; i++) {
			if ((i == N - 2))
				break;

			for (j = 2; j < N - 1; j++) {
				if (i == 0)
					break;

				/*Calculate r*/
				if (iterations < 10)
					phi_generate(0);
				else
					phi_generate(scheme);
				//Calculate numerical solution for the following iteration
				temp = m_grid[i][j];
				m_grid[i][j] = (m_grid[i][j + 1] * A_e + m_grid[i][j - 1] * A_w + m_grid[i][j - 2] * A_ww + m_grid[i - 1][j] * A_n + m_grid[i + 1][j] * A_s + m_grid[i + 2][j] * A_ss) / A_p;
				temp2 = m_grid[i][j];
				m_grid[i][j] = temp + SOR * (m_grid[i][j] - temp);
				residual = residual + abs(m_grid[i][j] - temp);
			}
		}

		MPI_Allreduce(&residual, &residual_after, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		residual = residual_after;

		if (iterations % BC_exchange_rate == 0) {
			for (i = 0; i < size; i++) {
				MPI_Bcast(&m_grid[partition[i]][0], (partition[rank + 1] - partition[rank]) * N, MPI_DOUBLE, i, MPI_COMM_WORLD);
			}
		}
	}
	time2 = MPI_Wtime() - time1;

	if (rank == 0)
		printinfo();
	/*-------------*/
	MPI_Finalize();
	return 0;
}

void printinfo() {
	write(&m_grid[0][0], N, N);
	cout << "-----------" << endl;
	cout << "Peclet number:	" << Pe << endl;
	cout << "Mesh number:" << N << "		Scheme:" << scheme << endl;
	cout << "SOR:	" << SOR << "		BC_exchange_rate:" << BC_exchange_rate << endl;
	cout << "Iterations:" << iterations << "		residual:" << residual << endl;
	cout << "Time:" << time2 << endl;
	cout << "-----------" << endl;
	return;
}

void phi_generate(int scheme_implement) {
	/*Calculate r*/
	r_e = (m_grid[i][j + 1] - m_grid[i][j]) / (m_grid[i][j] - m_grid[i][j - 1]);
	r_n = (m_grid[i - 1][j] - m_grid[i][j]) / (m_grid[i][j] - m_grid[i + 1][j]);
	r_s = (m_grid[i][j] - m_grid[i + 1][j]) / (m_grid[i + 1][j] - m_grid[i + 2][j]);
	r_w = (m_grid[i][j] - m_grid[i][j - 1]) / (m_grid[i][j - 1] - m_grid[i][j - 2]);

	/*Calculate phi*/
	switch (scheme_implement) {
	case 0:
		phi_e = 0, phi_n = 0, phi_s = 0, phi_w = 0;
		break;
	case 1:
		phi_e = 1, phi_n = 1, phi_s = 1, phi_w = 1;
		break;
	case 2:
		phi_e = r_e, phi_n = r_n, phi_s = r_s, phi_w = r_w;
		break;
	case 3:
		phi_e = max(0.0, min(2 * r_e, min(2.0, 2 * (r_e + abs(r_e)) / (r_e + 3))));
		phi_n = max(0.0, min(2 * r_n, min(2.0, 2 * (r_n + abs(r_n)) / (r_n + 3))));
		phi_s = max(0.0, min(2 * r_s, min(2.0, 2 * (r_s + abs(r_s)) / (r_s + 3))));
		phi_w = max(0.0, min(2 * r_w, min(2.0, 2 * (r_w + abs(r_w)) / (r_w + 3))));
		break;
	case 4:
		phi_e = max(0.0, min(r_e, 1.0));
		phi_n = max(0.0, min(r_n, 1.0));
		phi_s = max(0.0, min(r_s, 1.0));
		phi_w = max(0.0, min(r_w, 1.0));
		break;

	}

	A_e = D_e;
	A_w = ((m_rate * phi_e + m_rate * phi_w) * 0.5 + m_rate) + D_w;
	A_n = D_n;
	A_s = ((m_rate * phi_n + m_rate * phi_s) * 0.5 + m_rate) + D_s;
	A_ww = -m_rate * phi_w * 0.5;
	A_ss = -m_rate * phi_s * 0.5;
	A_p = A_e + A_w + A_n + A_s + A_ww + A_ss;
	return;
}

void write(double* a, int x, int y) {
	ofstream myfile("result.csv");
	int i, j;
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
			myfile << *(a + i * x + j) << ",";
		}
		myfile << endl;
	}

	myfile.close();
	return;
}
