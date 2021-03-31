#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI 3.14159265358979323846

using namespace std;
int N_el = 0, N_knot = 0;
FILE* in;

struct knot
{
	double x;
};

struct elem
{
	vector<knot> knots;
	vector<int> vertex_glob;
	int n_lambda, n_f;
	double gamma = 0.;
};

struct localA
{
	vector<vector<double>> locA;
};

struct slae
{
	vector<double> di;
	vector<double> al;
	vector<double> au;
	vector<int> ia;
	vector<double> b;
	vector<double> q;
};

double Lambda(int n_lambda, knot& kn)
{
	switch (n_lambda)
	{
	case(1):
		return 1.0;
	case(2):
		return exp(kn.x);
	}
}

double Sigma(int n_lambda, knot& kn)
{
	switch (n_lambda)
	{
	case(1):
		return 1.0;
	case(2):
		return kn.x;
	}
}

double Func(int n_f, knot& kn)
{
	switch (n_f)
	{
	case(1):
		return -2;
	case(2):
		return PI * exp(kn.x) * (PI * sin(PI * kn.x) - cos(PI * kn.x));
	}
}

double S1(int n_s)
{
	return 0.0;
}

double S2(int n_s)
{
	switch (n_s)
	{
	case(1):
		return 2;
	case(2):
		return -PI * exp(1);
	}
}

void Input_info(vector<double>& timeMesh, double a, double b)
{
	if (fopen_s(&in, "info.txt", "r") == 0)
	{
		int n = 0;
		double a, b;
		fscanf_s(in, "%d %lf %lf %d", &N_el, &a, &b, &n);
		timeMesh.resize(n);
		for (int i = 0; i < n; i++) {
			fscanf_s(in, "%fl", &timeMesh[i]);
		}
	}
	fclose(in);
}

void GetMesh(double& a, double& b)
{
	ofstream x("x.txt"), elems("elem.txt");
	double h = (b - a) / N_el, p = 0;
	x << N_el + 1 << endl;
	for (int i = 0; i <= N_el; i++)
	{
		p = a + h * i;
		x << p << endl;
		elems << i + 1 << " " << i + 2 << endl;
	}
	x.close();
	elems.close();
}

void Input_x(vector<knot>& knots)
{
	if (fopen_s(&in, "x.txt", "r") == 0)
	{
		fscanf_s(in, "%d", &N_knot);
		knots.resize(N_knot);
		for (int i = 0; i < N_knot; i++)
			fscanf_s(in, "%lf", &knots[i].x);
	}
	fclose(in);
}

void Input_elem(vector<elem>& elems, vector<knot>& knots)
{
	elems.resize(N_el);
	if (fopen_s(&in, "elem.txt", "r") == 0)
	{
		for (int i = 0; i < N_el; i++)
		{
			elems[i].knots.resize(2);
			elems[i].vertex_glob.resize(2);

			for (int j = 0; j < 2; j++)
			{
				fscanf_s(in, "%d", &elems[i].vertex_glob[j]);
				elems[i].vertex_glob[j]--;
				elems[i].knots[j] = knots[elems[i].vertex_glob[j]];
			}
		}
	}
	fclose(in);
}

void InputData(vector<double>& timeMesh, vector<knot>& knots, vector<elem>& elems) {
	double a = 0, b = 0;
	Input_info(timeMesh, a, b);

	GetMesh(a, b);
	Input_x(knots);
	Input_elem(elems, knots);
}

void CreatePortrait(slae A)
{
	A.di.resize(N_knot);
	A.ia.resize(N_knot);
	A.al.resize(N_knot - 1);
	A.au.resize(N_knot - 1);
	A.q.resize(N_knot);
	A.b.resize(N_knot);

	for (int i = 0; i < N_knot - 1; i++)
	{
		A.ia[i + 1] = i;
	}
}

void CreateLocal_b(vector<vector<double>>& bloc, vector<elem>& elems)
{
	double h = 0., f1 = 0., f2 = 0.;
	for (int i = 0; i < N_el; i++)
	{
		bloc[i].resize(2);

		f1 = Func(1, elems[i].knots[0]);
		f2 = Func(1, elems[i].knots[1]);

		h = elems[i].knots[1].x - elems[i].knots[0].x;
		bloc[i][0] = h / 6 * (2 * f1 + f2);
		bloc[i][1] = h / 6 * (f1 + 2 * f2);
	}
}

void Global_b(vector<double>& glob_b, vector<elem>& elems, vector<vector<double>>& bloc)
{
	for (int k = 0; k < N_el; k++)
		for (int iloc = 0; iloc < 2; iloc++)
		{
			int iglob = elems[k].vertex_glob[iloc];
			glob_b[iglob] += bloc[k][iloc];
		}
}

void CreateLocalA(vector<localA>& A, vector<elem>& elems)
{
	vector<vector<double>> G(2), M(2);
	double h = 0., l1 = 0., l2 = 0., sigma1 = 0., sigma2 = 0.;
	for (int i = 0; i < 2; i++)
	{
		G[i].resize(2);
		M[i].resize(2);
	}
	for (int i = 0; i < N_el; i++)
	{
		A[i].locA.resize(2);
		for (int j = 0; j < 2; j++)
		{
			A[i].locA[j].resize(2);
		}

		l1 = Lambda(1, elems[i].knots[0]);
		l2 = Lambda(1, elems[i].knots[1]);

		sigma1 = Sigma(1, elems[i].knots[0]);
		sigma2 = Sigma(1, elems[i].knots[1]);

		h = elems[i].knots[1].x - elems[i].knots[0].x;
		G[0][0] = G[1][1] = (l1 + l2) / (2 * h);
		G[0][1] = G[1][0] = -(l1 + l2) / (2 * h);

		M[0][0] = M[1][1] = (sigma1 + sigma2) * h / 3;
		M[0][1] = M[1][0] = (sigma1 + sigma2) * h / 6;

		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				A[i].locA[j][k] = M[j][k] + G[j][k];
	}
	for (int i = 0; i < 2; i++)
	{
		double sum = 0.;
		for (int j = 0; j < 2; j++)
		{
			sum += G[i][j];
		}
		cout << sum << endl;
	}
}

void GlobalA(vector<vector<double>>& globA, vector<localA>& A, vector<elem>& elems)
{

}

void GlobPlusCond1(vector<vector<double>>& globA, vector<double>& glob_b)
{
	globA[0][0] = 1;
	for (int j = 1; j < N_knot; j++)
	{
		globA[0][j] = 0;
	}
	glob_b[0] = S1(1);
}

void GlobPlusCond2(vector<double>& glob_b)
{
	glob_b[N_knot - 1] += S2(1);
}

void Solution(vector<double>& di, vector<double>& al, vector<double>& au, vector<double>& F, vector<double>& q)
{
	q.resize(N_knot);
	double buf = 0.;
	vector<double> alpha(N_knot - 1), beta(N_knot);
	alpha[0] = -au[0] / di[0];
	beta[0] = F[0] / di[0];

	for (int i = 1; i < N_knot - 1; i++)
	{
		buf = di[i] + al[i - 1] * alpha[i - 1];
		alpha[i] = -au[i] / buf;
		beta[i] = (F[i] - al[i - 1] * beta[i - 1]) / buf;
	}
	q[N_knot - 1] = (F[N_knot - 1] - al[N_knot - 2] * beta[N_knot - 2]) / (di[N_knot - 1] + al[N_knot - 2] * alpha[N_knot - 2]);
	for (int i = N_knot - 2; i >= 0; i--)
	{
		q[i] = alpha[i] * q[i + 1] + beta[i];
	}
}

void Gauss(vector < vector < double>>& A, vector<double>& res, vector<double>& F)
{
	int n = res.size(), i = 0;
	double max = 0, buf = 0;

	for (int k = 0; k < n; k++)
	{
		max = abs(A[k][k]);
		i = k;
		for (int m = k + 1; m < n; m++)
		{
			buf = A[m][k];
			if (abs(buf) > max)
			{
				i = m;
				max = abs(buf);
			}
		}

		if (max == 0)
			cout << "System havent solution" << endl;

		if (i != k)
		{
			A[i].swap(A[k]);
			buf = F[i];
			F[i] = F[k];
			F[k] = buf;
		}

		max = A[k][k];
		A[k][k] = 1;

		for (int j = k + 1; j < n; j++)
			A[k][j] /= max;

		F[k] /= max;

		for (i = k + 1; i < n; i++)
		{
			buf = A[i][k];
			A[i][k] = 0;
			if (buf != 0)
			{
				for (int j = k + 1; j < n; j++)
					A[i][j] -= buf * A[k][j];
				F[i] -= buf * F[k];
			}
		}
	}
	for (i = n - 1; i >= 0; i--)
	{
		buf = 0;
		for (int j = n - 1; j > i; j--)
			buf += A[i][j] * res[j];
		res[i] = F[i] - buf;
	}
}

double Norm(vector<double> x)
{
	double res = 0.;
	int n = x.size();
	for (int i = 0; i < n; i++)
		res += x[i] * x[i];
	return sqrt(res);
}

int main()
{
	vector<double> timeMesh;
	vector<knot> knots;
	vector<elem> elems;

	InputData(timeMesh, knots, elems);

	slae A;
	CreatePortrait(A);
	vector<localA> Aloc(N_el);
	vector<vector<double>> bloc(N_el);
	vector<double> glob_b(N_knot);

	CreateLocal_b(bloc, elems);
	Global_b(glob_b, elems, bloc);
	CreateLocalA(Aloc, elems);

	//Gauss(globA, q, glob_b);

	//FILE* out;
	//fopen_s(&out, "out.txt", "w");
	//fprintf_s(out, "%7s", "x");
	//fprintf_s(out, "%2s", " |");
	//fprintf_s(out, "%15s", "Yn*");
	//fprintf_s(out, "%20s", "Yn\n");
	//for (int i = 0; i < 55; i++)
	//	fprintf_s(out, "%1s", "-");
	//fprintf_s(out, "%s", "\n");

	//vector<double> buf(N_knot);
	//for (int i = 0; i < N_knot; i++)
	//{
	//	buf[i] = sin(knots[i].x*PI);
	//	fprintf_s(out, "%7f", knots[i].x);
	//	fprintf_s(out, "%2s", " |");
	//	fprintf_s(out, "%.14e ", buf[i]);
	//	fprintf_s(out, " %.14e ", q[i]);
	//	buf[i] -= q[i];
	//	fprintf_s(out, "%s", "\n");
	//}
	//fprintf_s(out, "%s", "\n");
	//fprintf_s(out, "%5s", "||Yn*-Yn|| = ");
	//fprintf_s(out, "%.14e ", Norm(buf));
	//fclose(out);

	return 0;
}