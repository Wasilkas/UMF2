#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define PI 3.14159265358979323846

using namespace std;
int N_el = 0, N_knot = 0;
double EPSILON = 1e-13, DELTA = 1e-13;
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

struct slae
{
	vector<double> di;
	vector<double> al;
	vector<double> au;
	vector<int> ia;
	vector<double> b;
	vector<vector<double>> q;
};

struct matrix
{
	vector<double> di;
	vector<double> al;
	vector<double> au;
	vector<int> ia;

	matrix() {
		di.resize(N_knot);
		ia.resize(N_knot + 1);

		ia[0] = 0;
		for (int i = 0; i < N_knot; i++)
			ia[i + 1] = i;

		au.resize(N_knot - 1);
		al.resize(N_knot - 1);
	}
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

double Sigma(int n_lambda, knot& kn, double q, double t)
{
	switch (n_lambda)
	{
	case(1):
		return kn.x*t;
	case(2):
		return kn.x;
	}
}

double Func(int n_f, knot& kn, double t)
{
	switch (n_f)
	{
	case(1):
		return 5 * kn.x + 3;
	case(2):
		return 2 * t * kn.x;
	case(3):
		return kn.x * kn.x * t;
	}
}

double Q0(int n_f, knot& kn)
{
	switch (n_f)
	{
	case(1):
		return 0;
	}
}

double S1(int n_s, double t)
{
	switch (n_s)
	{
	case(1):
		return 3 * t;
	case(2):
		return 53 * t;
	case(3):
		return 0;
	case(4):
		return 10 * t;
	}
}

double S2(int n_s, double t)
{
	switch (n_s)
	{
	case(1):
		return 50 * t + 30;
	case(2):
		return -PI * exp(1);
	}
}

void Input_info(vector<double>& timeMesh, double &a, double &b)
{
	int n = 0;
	if (fopen_s(&in, "info.txt", "r") == 0)
	{
		fscanf_s(in, "%d %lf %lf %d", &N_el, &a, &b, &n);
		timeMesh.resize(n);
		for (int i = 0; i < n; i++) {
			fscanf_s(in, "%lf", &timeMesh[i]);
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

void CreatePortrait(slae &A, int timeSize)
{
	A.di.resize(N_knot);
	A.ia.resize(N_knot + 1);
	A.al.resize(N_knot - 1);
	A.au.resize(N_knot - 1);
	A.q.resize(timeSize);

	for (int i = 0; i < timeSize; i++)
	{
		A.q[i].resize(N_knot);
	}

	A.b.resize(N_knot);

	for (int i = 0; i < N_knot; i++)
	{
		A.ia[i + 1] = i;
	}
}

void Mult(vector<int>& ia, vector<double>& al, vector<double>& au, vector<double>& di, vector<double>& vec, vector<double> &res) {
	int dim = vec.size(), j = 0;
	res.resize(dim);
	for (int i = 0; i < dim; i++) {
		res[i] = di[i] * vec[i];
		j = i - ia[i + 1] + ia[i];
		for (int k = ia[i]; k < ia[i + 1]; j++, k++) {
			res[i] += al[k] * vec[j];
			res[j] += au[k] * vec[i];
		}
	}
}

void CreateGlobalM(matrix& MG, vector<elem>& elems, double &t, vector<double> &q)
{
	vector<vector<double>> M(2);
	double h = 0., sigma1 = 0., sigma2 = 0;
	for (int i = 0; i < 2; i++)
	{
		M[i].resize(2);
	}
	for (int i = 0; i < N_el; i++)
	{
		sigma1 = Sigma(1, elems[i].knots[0], q[elems[i].vertex_glob[0]], t);
		sigma2 = Sigma(1, elems[i].knots[1], q[elems[i].vertex_glob[1]], t);

		h = elems[i].knots[1].x - elems[i].knots[0].x;

		M[0][0] = (3*sigma1 + sigma2) * h / 12;
		M[1][1] = (sigma1 + 3 * sigma2) * h / 12;
		M[0][1] = M[1][0] = (sigma1 + sigma2)/2 * h / 6;

		MG.di[elems[i].vertex_glob[0]] += M[0][0];
		MG.di[elems[i].vertex_glob[1]] += M[1][1];
		MG.au[i] += M[0][1];
		MG.al[i] += M[1][0];
	}
}

void CreateGlobalG(matrix& GG, vector<elem>& elems, double& t) {
	vector<vector<double>> G(2);
	double h = 0., l1 = 0., l2 = 0.;

	G[0].resize(2);
	G[1].resize(2);

	for (int i = 0; i < N_el; i++)
	{
		l1 = Lambda(1, elems[i].knots[0]);
		l2 = Lambda(1, elems[i].knots[1]);

		h = elems[i].knots[1].x - elems[i].knots[0].x;
		G[0][0] = G[1][1] = (l1 + l2) / (2 * h);
		G[0][1] = G[1][0] = -(l1 + l2) / (2 * h);

		GG.di[elems[i].vertex_glob[0]] += G[0][0];
		GG.di[elems[i].vertex_glob[1]] += G[1][1];
		GG.au[i] += G[0][1];
		GG.al[i] += G[1][0];
	}


}

void CreateGlobalB(vector<double> &bg, vector<elem>& elems, double& t) {
	double h = 0., f1 = 0., f2 = 0.;
	vector<double> b(2);

	for (int i = 0; i < N_el; i++)
	{
		f1 = Func(3, elems[i].knots[0], t);
		f2 = Func(3, elems[i].knots[1], t);

		h = elems[i].knots[1].x - elems[i].knots[0].x;

		b[0] = (h / 6) * (2 * f1 + f2);
		b[1] = (h / 6) * (f1 + 2 * f2);

		bg[elems[i].vertex_glob[0]] += b[0];
		bg[elems[i].vertex_glob[1]] += b[1];
	}
}

void CreateGlobalA(slae& A, matrix& M, matrix& G, vector<double>& b, vector<double> &q0, double &deltaT) {
	vector<double> Mq(N_knot);
	Mult(M.ia, M.al, M.au, M.di, q0, Mq);
	
	for (int i = 0; i < N_knot; i++) {
		A.di[i] = M.di[i] / deltaT + G.di[i];
		A.b[i] = b[i] + Mq[i] / deltaT;
	}

	for (int i = 0; i < N_knot - 1; i++) {
		A.al[i] = M.al[i] / deltaT + G.al[i];
		A.au[i] = M.au[i] / deltaT + G.au[i];
	}
}

void calcQ0(vector<double> &q0, vector<knot> &knots) 
{
	for (int i = 0; i < N_knot; i++)
	{
		q0[i] = Q0(1, knots[i]);
	}
}

void GlobPlusCond1(slae& A, double &t)
{
	A.di[0] = 1;
	A.au[0] = 0;
	A.b[0] = S1(3, t);
	A.di[N_knot - 1] = 1;
	A.al[A.al.size() - 1] = 0;
	A.b[N_knot - 1] = S1(4, t);
}

void GlobPlusCond2(vector<double>& b, double &t)
{
	b[N_knot - 1] += S2(1, t);
}

double Norm(vector<double> x)
{
	double res = 0.;
	int n = x.size();
	for (int i = 0; i < n; i++)
		res += x[i] * x[i];
	return sqrt(res);
}

void LUwithStar(vector<int>& ia, vector<double>& al, vector<double>& au, vector<double>& di) {
	int dim = di.size(), j = 0, ik0 = 0, i_st = 0, j_st = 0, jk0 = 0;
	double sumL = 0., sumU = 0., sumDi = 0.;
	for (int i = 0; i < dim; i++) {
		sumDi = 0.;
		j = i - (ia[i + 1] - ia[i]);
		i_st = j;
			for (int k = ia[i]; k < ia[i + 1]; k++, j++) {
				sumL = 0.;
				sumU = 0.;
				j_st = j - (ia[j + 1] - ia[j]);
				if (i_st > j_st) {
					ik0 = ia[i];
					jk0 = ia[j] + (i_st - j_st);
				}
				else {
					ik0 = ia[i] + (j_st - i_st);
					jk0 = ia[j];
				}
				for (int ik = ik0; ik < k; ik++, jk0++) {
					sumL += al[ik] * au[jk0];
					sumU += au[ik] * al[jk0];
				}
				al[k] = (al[k] - sumL) / di[j];
				au[k] -= sumU;
				sumDi += al[k] * au[k];
			}
		di[i] -= sumDi;
	}
}

void ForwardMotion(vector<int>& ia, vector<double>& al, vector<double>& vec) {
	int dim = vec.size(), j = 0, i0 = 0, i1 = 0;
	double sum = 0.;
	for (int i = 1; i < dim; i++) {
		sum = 0.;
		i0 = ia[i];
		i1 = ia[i + 1];
		j = i - (i1 - i0);
		for (int k = i0; k < i1; k++, j++) {
			sum += al[k] * vec[j];
		}
		vec[i] = vec[i] - sum;
	}
}
void ReverseMotion(vector<int>& ia, vector<double>& au, vector<double>& di, vector<double>& vec) {
	int dim = vec.size(), j = 0;
	double xi = 0;
	for (int i = dim - 1; i >= 0; i--) {
		xi = vec[i] / di[i];
		j = i - 1;
		for (int k = ia[i + 1] - 1; k >= ia[i]; k--, j--) {
			vec[j] -= au[k] * xi;
		}
		vec[i] = xi;
	}
}

double CalcResidual(slae& A, vector<double> &q, vector<double> &b)
{
	vector<double> v(q.size());
	Mult(A.ia, A.al, A.au, A.di, q, v);

	for (int i = 0; i < v.size(); i++)
	{
		v[i] -= b[i];
	}

	return Norm(v) / Norm(b);
}

double CalcChange(vector<double>& q_prev, const vector<double>& q)
{
	for (int i = 0; i < q.size(); i++)
		q_prev[i] = q[i] - q_prev[i];

	return Norm(q_prev) / Norm(q);
}

int main()
{
	vector<double> timeMesh;
	vector<knot> knots;
	vector<elem> elems;

	InputData(timeMesh, knots, elems);

	slae A;
	matrix M, G;
	vector<double> b(N_knot), b_q0(N_knot), q_prev(N_knot);
	
	CreatePortrait(A, timeMesh.size());

	calcQ0(A.q[0], knots);

	cout << "q" << 0 << "C" << ":\t";
	for (int i = 0; i < A.q[0].size(); i++)
	{
		cout << A.q[0][i] << "\t";
	}
	cout << endl;
	cout << "q" << 0 << "A" << ":\t";
	for (int i = 0; i < A.q[0].size(); i++)
	{
		cout << 0 << "\t";
	}
	cout << endl;

	double deltaT = 0., t = 0., resid = EPSILON + 1, changeQ = DELTA + 1;
	for (int i = 1; i < timeMesh.size(); i++)
	{
		t = timeMesh[i];
		//q_prev[0] = 0;
		//q_prev[1] = 5 * t;
		//q_prev[2] = 10 * t;
		deltaT = t - timeMesh[i - 1];
		CreateGlobalM(M, elems, t, q_prev);
		CreateGlobalG(G, elems, t);
		CreateGlobalB(b, elems, t);
		CreateGlobalA(A, M, G, b, A.q[i - 1], deltaT);
		GlobPlusCond1(A, t);
		/*b_q0 = A.b;
		resid = CalcResidual(A, q_prev, b_q0);

		for (int iter = 1; iter < 1000 && resid > EPSILON && changeQ > DELTA; iter++)*/
		//{
			LUwithStar(A.ia, A.al, A.au, A.di);
			ForwardMotion(A.ia, A.al, A.b);
			ReverseMotion(A.ia, A.au, A.di, A.b);
		//	
		//	changeQ = CalcChange(q_prev, A.b);
		//	q_prev = A.b;

		//	CreateGlobalM(M, elems, t, q_prev);
		//	CreateGlobalA(A, M, G, b, A.q[i - 1], deltaT);
		//	GlobPlusCond1(A, t);

		//	resid = CalcResidual(A, q_prev, b_q0);

		//	cout << iter << "\t\t" << "Residual: " << resid << "\t\t" << "changeQ: " << changeQ << endl;
		//}
		//cout << endl;

		A.q[i] = A.b;
		cout << "q" << i << "C" << ":\t";
		for (int j = 0; j < A.q[0].size(); j++)
		{
			cout << A.q[i][j] << "\t";
		}
		cout << endl;
		cout << "q" << i << "A" << ":\t";
		for (int j = 0; j < A.q[0].size(); j++)
		{
			//cout << 5*knots[j].x*t + 3*t << "\t";
			cout << knots[j].x * t << "\t";
		}
		cout << "\n\n\n";
	}

	return 0;
}