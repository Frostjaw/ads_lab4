// lab 4_2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <ctime>
#include <clocale>
using namespace std;
#define MAX 1.7976931348623158e+308

struct coordinate {
	int x,y;
};

class wfgraph {
	int n;
	double** mat;
	int** smejmat;  //������� ��������� ��� ��� ������
	int* obhod;

	//��� ������ � ������
	int* optsol;
	int* check;
	int ncur;
public:
	double wminost;
	double wobhod;
	double weight_optsol;
	wfgraph();
	wfgraph(int n2);
	~wfgraph();
	void get_span_tree();
	double distance(coordinate a, coordinate b);
	void input_edges(coordinate* t);
	void add_edge(int a, int b);
	//void randgraph();
	void deep(int curver, int *R, int *r2, int &curnum);
	int* DFS();
	void komminostov();
	void branch_and_bound();
	double weight_obhoda(int* obhod, int n);
	void first_sol();
	void find_optsol(int* s);
	//friend ostream& operator<<(ostream &stream, wfgraph &t);
	void coutmat();
	void coutsmejmat();
	void coutobhod();
	void cout_optsol();
};

wfgraph::wfgraph() {
	n = 0; 
	mat = NULL; smejmat = NULL; obhod = NULL; optsol = NULL; check = NULL;
}

wfgraph::wfgraph(int n2) {
	int i, j; 
	n = n2;
	mat = new double*[n]; smejmat = new int*[n]; obhod = new int[n+1];
	for (i = 0;i < n;i++) { mat[i] = new double[n]; smejmat[i] = new int[n]; }
	for (i = 0;i < n;i++)
		for (j = 0;j < n;j++) smejmat[i][j] = 0;
}

wfgraph::~wfgraph() {
	for (int i = 0;i < n;i++) {
		delete[] mat[i]; delete[] smejmat[i];
	}
	delete[] mat; delete[] smejmat; delete[] obhod; delete[] optsol; delete[] check;
}

void wfgraph::get_span_tree() {
	double wmin;
	wminost = 0; //������������� wminost ��� ���������� � add_edge
	int i, j, vm, *b = new int[n];
	b[0] = -1;
	for (i = 1; i < n; i++) b[i] = 0;
	for (i = 1; i < n; i++) {
		wmin = MAX; vm = 0;
		for (j = 1; j < n; j++)
			if (b[j] != -1 && wmin > mat[j][b[j]])
			{
				vm = j; wmin = mat[j][b[j]];
			}
		if (!vm) return;
		add_edge(vm, b[vm]); b[vm] = -1;
		for (j = 1; j < n; j++)
			if (b[j] != -1 && mat[j][b[j]] > mat[j][vm]) b[j] = vm;
	}
}

//void wfgraph::randgraph() {
//	for (int i = 0;i < n;i++)
//		for (int j = 0;j < n;j++)
//			mat[i][j] = rand() % +20;
//	for (int i = 0;i < n;i++) mat[i][i] = MAX;
//}

double wfgraph::distance(coordinate a, coordinate b) {
	return sqrt((a.x-b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

void wfgraph::input_edges(coordinate* t) {
	for (int i = 0;i < n;i++)
		for (int j = 0;j < n;j++) { 
			if (i == j) mat[i][j] = MAX;
			else mat[i][j] = distance(t[i], t[j]);
		}
}

void wfgraph::add_edge(int a, int b) { //+ ���������� wminost
	smejmat[a][b] = 1; smejmat[b][a] = 1; wminost += mat[a][b];
}

//ostream& operator<<(ostream &stream, wfgraph &t) {
//	for (int i = 0;i < t.n;i++) {
//		for (int j = 0;j < t.n;j++) {
//			stream.width(3); stream.precision(3);
//			if (t.mat[i][j] == MAX) stream << "inf" << ' ';
//			else stream << round(t.mat[i][j]) << ' ';
//		}
//		stream << endl;
//	}
//	return stream;
//}

void wfgraph::deep(int cver, int *R, int *r2, int &cnum) // ������ r2 - ������� ������� ��� ������ � �������
{
	r2[cnum] = cver; R[cver] = ++cnum;
	for (int i = 0; i < n; i++)
		if (smejmat[cver][i] && !R[i]) deep(i, R, r2, cnum);
}

int* wfgraph::DFS()
{
	int i, cnum, *R = new int[n], *r2 = new int[n+1];
	for (i = 0; i < n; i++) R[i] = 0;
	for (cnum = i = 0; i < n; i++)
		if (!R[i]) deep(i, R, r2, cnum);
	return r2;
}

double wfgraph::weight_obhoda(int* obhod, int n) { // n - ���������� ������ � ������ (+1) 
	double temp = 0;
	for (int i = 0;i < n-1;i++) temp += mat[obhod[i]][obhod[i + 1]];
	return temp;
}

//void wfgraph::deep(int cver, int *R, int &cnum)
//{
//	R[cver] = ++cnum;
//	for (int i = 0; i < n; i++)
//		if (smejmat[cver][i] && !R[i]) deep(i, R, cnum);
//}

//int* wfgraph::DFS()
//{
//	int i, cnum, *R = new int[n];
//	for (i = 0; i < n; i++) R[i] = 0;
//	for (cnum = i = 0; i < n; i++)
//		if (!R[i]) deep(i, R, cnum);
//	return R;
//}

void wfgraph::komminostov() {
	obhod = DFS(); obhod[n] = 0; // ����������� �������
	wobhod = weight_obhoda(obhod,n+1);
}

void wfgraph::first_sol() {
	int i;
	optsol = new int[n+1]; weight_optsol = 0;
	for (i = 0;i < n;i++) optsol[i] = i; optsol[n] = 0;
	for (i = 0;i < n;i++) weight_optsol += mat[optsol[i]][optsol[i + 1]];
}

void wfgraph::branch_and_bound() {
	int* sol = new int[n + 1]; // sol - "���������" �������, ��� ������� ���������
	check = new int[n];
	sol[0] = 0; check[0] = 1; ncur = 1;
	find_optsol(sol);
	delete[] sol;
}

void wfgraph::find_optsol(int* s) {
	double weight;
	for (int i = 1; i < n; i++) {
		if (check[i] != 1) {
			s[ncur] = i; s[ncur + 1] = 0; check[i] = 1; ncur++; //��������� ������� i � �������� ������� + ����������� ������� �������
			weight = weight_obhoda(s, ncur+1); 
			if (ncur == n && weight < weight_optsol) //�������� �� ������ ������� � weight < weight_optsol
			{
				for (int i = 0; i <= n; i++) optsol[i] = s[i];
				weight_optsol = weight;
			}
			else if (weight < weight_optsol) find_optsol(s);
			ncur--; check[i] = 0; //������� ������� i �� �������� �������
		}
	}
}

void wfgraph::coutmat() {
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout.width(3); cout.precision(3);
			if (mat[i][j] == MAX) cout << "inf" << ' ';
			else cout << round(mat[i][j]) << ' ';
		}
		cout << endl;
	}
}

void wfgraph::coutsmejmat() {
	for (int i = 0;i < n;i++) {
		for (int j = 0;j < n;j++) {
			cout.width(3); cout.precision(3);
			if (mat[i][j] == MAX) cout << "inf" << ' ';
			else cout << round(smejmat[i][j]) << ' ';
		}
		cout << endl;
	}
}

void wfgraph::coutobhod() {
	for (int i = 0;i <= n;i++) cout << obhod[i] + 1 << ' ';           // + 1 ��� ��������� ������
}

void wfgraph::cout_optsol() {
	for (int i=0;i<=n;i++) cout << optsol[i] +1 << ' ';
}

int main()
{
	setlocale(LC_ALL, "Russian");
	srand(time(0));
	/*wfgraph test(5);
	test.randgraph();
	cout << test;*/

	//���� ������� ���������
	int nver,i; double t1, t2;

	/*cout << "������� ���-�� ������" << endl; cin >> nver;*/
	cout << "��� ��������� ��� 10, 15 ������: " << endl;
	// 10, 15 ������ ��� ����� ����������
	for (int i = 10; i <= 15; i += 5) {
		nver = i; cout << "������: " << nver << endl;
		coordinate* testmap;
		testmap = new coordinate[nver];
		for (i = 0; i < nver; i++) { testmap[i].x = rand() % +15; testmap[i].y = rand() % +15; }
		cout << "�������� �����: " << endl;
		for (i = 0; i < nver; i++) cout << i + 1 << ":<" << testmap[i].x << ',' << testmap[i].y << '>' << ' ';
		cout << endl << endl;
		wfgraph test(nver);
		test.input_edges(testmap);

		//����������� �� ������ min ������
		t1 = clock();
		test.get_span_tree();
		test.komminostov();
		t2 = clock();
		cout << "������� �����:" << endl;
		test.coutmat(); cout << endl;
		cout << "������������ �������� ������ �� ������ ������������ ������: " << endl << endl;
		cout << "������� ��������� ������������ ������" << endl;
		test.coutsmejmat(); cout << endl;
		cout << "����������� �������: " << endl; test.coutobhod();
		cout << endl << "��� ������������ ������ = " << test.wminost << ' ' << endl;
		cout << "��� �������� �� ������ min ������ = " << test.wobhod << ' ' << endl;
		cout << "������ �������� = " << test.wobhod / test.wminost << endl;
		cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl;

		//���������� �������� � ����������� (����� ������ � ������)

		//����������� ������� (������)
		test.first_sol();
		cout << endl << "���������� �������� � �����������: " << endl;
		cout << endl << "������ ����������� �������: " << endl;
		test.cout_optsol(); cout << endl;
		cout << "��� = " << test.weight_optsol << endl;

		//����� ������ � ������
		t1 = clock();
		test.branch_and_bound();
		t2 = clock();
		cout << endl << "�������:" << endl;
		test.cout_optsol(); cout << endl;
		cout << "��� = " << test.weight_optsol << endl;
		cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl << endl;
	}

	// 50, 100, 500, 1000 ������ ��� �������������

	cout << "���������� �������� ��� 50, 100, 500, 1000 ������: " << endl << endl;
	// 50:
	nver = 50; cout << "������: " << nver << endl;
	coordinate* testmap50;
	testmap50 = new coordinate[nver];
	for (i = 0; i < nver; i++) { testmap50[i].x = rand() % +15; testmap50[i].y = rand() % +15; }
	wfgraph test50(nver);
	test50.input_edges(testmap50);
	//����������� �� ������ min ������
	t1 = clock();
	test50.get_span_tree();
	test50.komminostov();
	t2 = clock();
	cout << "��� ������������ ������ = " << test50.wminost << ' ' << endl;
	cout << "��� �������� �� ������ min ������ = " << test50.wobhod << ' ' << endl;
	cout << "������ �������� = " << test50.wobhod / test50.wminost << endl;
	cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl << endl;

	//100:
	nver = 100; cout << "������: " << nver << endl;
	coordinate* testmap100;
	testmap100 = new coordinate[nver];
	for (i = 0; i < nver; i++) { testmap100[i].x = rand() % +15; testmap100[i].y = rand() % +15; }
	wfgraph test100(nver);
	test100.input_edges(testmap100);
	//����������� �� ������ min ������
	t1 = clock();
	test100.get_span_tree();
	test100.komminostov();
	t2 = clock();
	cout << "��� ������������ ������ = " << test100.wminost << ' ' << endl;
	cout << "��� �������� �� ������ min ������ = " << test100.wobhod << ' ' << endl;
	cout << "������ �������� = " << test100.wobhod / test100.wminost << endl;
	cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl << endl;

	// 500:
	nver = 500; cout << "������: " << nver << endl;
	coordinate* testmap500;
	testmap500 = new coordinate[nver];
	for (i = 0; i < nver; i++) { testmap500[i].x = rand() % +15; testmap500[i].y = rand() % +15; }
	wfgraph test500(nver);
	test500.input_edges(testmap500);
	//����������� �� ������ min ������
	t1 = clock();
	test500.get_span_tree();
	test500.komminostov();
	t2 = clock();
	cout << "��� ������������ ������ = " << test500.wminost << ' ' << endl;
	cout << "��� �������� �� ������ min ������ = " << test500.wobhod << ' ' << endl;
	cout << "������ �������� = " << test500.wobhod / test500.wminost << endl;
	cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl << endl;

	// 1000:
	nver = 1000; cout << "������: " << nver << endl;
	coordinate* testmap1000;
	testmap1000 = new coordinate[nver];
	for (i = 0; i < nver; i++) { testmap1000[i].x = rand() % +15; testmap1000[i].y = rand() % +15; }
	wfgraph test1000(nver);
	test1000.input_edges(testmap1000);
	//����������� �� ������ min ������
	t1 = clock();
	test1000.get_span_tree();
	test1000.komminostov();
	t2 = clock();
	cout << "��� ������������ ������ = " << test1000.wminost << ' ' << endl;
	cout << "��� �������� �� ������ min ������ = " << test1000.wobhod << ' ' << endl;
	cout << "������ �������� = " << test1000.wobhod / test1000.wminost << endl;
	cout << "����� ������ = " << (t2 - t1) / 1000 << ' ' << endl << endl;

    return 0;
}

