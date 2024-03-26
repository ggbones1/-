#include <string>
#include <sstream> 
#include <fstream> 
#include <iomanip>
#include <iostream> 
#include<cmath>
#include<cstdlib>
//#include <windows.h>
#include<vector>
#include "common.h"
using namespace std;
//
//*************************************************************************************************
struct GridData
{
	int Type;
	double Den, Ems, Pt, P0, Dx, Dy, Dd, Ddx, Ddy, Px, Py;
	double Vx, Vy, EVx, EVy, MVx, MVy, Fx, Fy;
};

double(*Dist)[DY][DQ];
double(*Temp)[DY][DQ];
GridData(*Grid)[DY];
double CoDen[512][3];
int CoNum = 0;
char MacroName[128];
double Mass, Density, Tau, ReTau, Viscosity, Radius, MaxSpeed;
double DenG, DenL, G, A, B, R, K, T, Tr, Tc, Rc, Pc, Kr, Fe, Fm;
int    NowStep, AllStep, ShowStep, SaveStep, BeginTime, StepTime, LastTime, MModel, FModel, No, TaskNum = 1;
double Fn[2000], Dn[2000];
int flag = 0;

ofstream TrackFile;
std::ofstream File;

//*************************************************************************************************
void Initialize();
void SetMultiphase();
void SetFlowField();
void ShowData();
void SaveData();
void NonidealForce();
void GlobalCollide();
void MacroCalculate();
//int  TimeInterval() { StepTime = GetTickCount() - LastTime;  LastTime = GetTickCount();  return StepTime; }
void SaveDen() {
	//if (No<TaskNum-1) {
	if (Tr < 1.0) {
		CoDen[CoNum][0] = T / Tc;
		CoDen[CoNum][1] = Grid[DX / 2][0].Den / Rc;
		CoDen[CoNum][2] = Grid[DX / 2][DY / 2].Den / Rc;
		CoNum++;

	}
	else {
		string FileName = Name(MModel);
		FileName = "../results/T=0.6-TrDen--" + FileName + ".txt";
		ofstream File(FileName);
		File << "Density   Tr" << endl;
		for (int i = 0; i < CoNum; ++i)      File << CoDen[i][1] << "   " << CoDen[i][0] << endl;
		File << "1   1" << endl;
		for (int i = CoNum - 1; i >= 0; --i)   File << CoDen[i][2] << "   " << CoDen[i][0] << endl;
		File.close();
	}
}

//*************************************************************************************************
int main(int argc, char *argv[])
{

	for (Tr = 0.6; Tr < 1.0; Tr += 0.1)
	{
		/*for (No = 0; No < TaskNum; ++No)
		{*/
		Initialize();
		SetMultiphase();
		SetFlowField();
		ShowData();

		for (NowStep = 1; NowStep <= AllStep; ++NowStep)
		{
			NonidealForce();
			GlobalCollide();
			MacroCalculate();

			if (NowStep%SaveStep == 0) { SaveData(); }
			if (NowStep%ShowStep == 0) {//TimeInterval(); 
				ShowData();
			}
		}
		SaveDen();
	}
	if (Dist != 0)	delete[] Dist;
	if (Temp != 0)	delete[] Temp;
	if (Grid != 0)	delete[] Grid;
	SaveDen();
	return 0;
}

//*************************************************************************************************
void Initialize()
{
	/*TaskNum = 12;
	Tr = 0.5;
	flag = 0;*/
	Mass = 0;
	Density = 1;
	G = -1;
	MModel = MP_P0PRW;
	FModel = FI_KUPER;
	Tau = 1.2;
	Kr = 0.1;

	AllStep = 2 * 10000;
	SaveStep = 1000;
	ShowStep = 1000;
	NowStep = StepTime = 0;
	//BeginTime = LastTime = GetTickCount();

	Dist = new double[DX][DY][DQ];
	Temp = new double[DX][DY][DQ];
	Grid = new GridData[DX][DY];

	if (Grid == 0 || Dist == 0 || Temp == 0)
	{
		cout << "Memory allocation is error ..." << endl << flush;
		cin.get();
		exit(1);
	}
}

//*************************************************************************************************
void SetFlowField()
{
	ReTau = double(1) / Tau;
	Viscosity = (Tau * 2 - 1) / 6;
	FOR_iDX_jDY
	{
		GridData & t = Grid[i][j];

		t.Type = FLUID;

		double j1 = D(j) - DY / 4, j2 = D(j) - DY / 4 * 3, Width = 5;
		t.Den = DenG + (DenL - DenG) / 2 * (tanh(j1 * 2 / Width) - tanh(j2 * 2 / Width));
		for (int f = 0; f < DQ; ++f)
		{
			Dist[i][j][f] = Temp[i][j][f] = Feq(f, t.Den, 0, 0);
		}
	}

	FOR_iDX_jDY_Fluid  Mass += Grid[i][j].Den;

	cout << endl << " Multiphase: " << Name(MModel) << "   DX=" << DX << "   DY=" << DY << "   Tau=" << setprecision(2) << Tau << "   Tr=" << Tr << "   DenG=" << setprecision(5) << DenG << "   DenL=" << DenL << endl;
	cout << "**************************************************************************************" << endl;
}
//******************************************************************************************************************************
void SetMultiphase()
{
	K = 0.001;
	switch (MModel)
	{
	case MP_P0VDW:
		A = D(9) / 49;   B = D(2) / 21;   Tc = D(4) / 7;     Rc = D(7) / 2;
		break;

	case MP_P0CSE:
		A = 1.0;       B = 1.0;   Tc = 0.09432870;     Rc = 0.13044388;
		break;

	case MP_P0RKE:
	case MP_P0RKS: {
		A = D(2) / 49;   B = D(2) / 21;   Rc = 2.729171;  Tc = (MModel == MP_P0RKE ? 0.196133 : 0.086861); double  w = 0.344;
		A *= (MModel == MP_P0RKE ? D(1) / sqrt(Tc*Tr) : Sq(D(1) + (0.480 + 1.574*w - 0.176*w*w)*(D(1) - sqrt(Tr)))); }
				   break;

	case MP_P0PRM:
	case MP_P0PRW: {
		A = D(2) / 49;   B = D(2) / 21;   Tc = 0.072919;   Rc = 2.65730416; double w = (MModel == MP_P0PRW) ? 0.344 : 0.011;
		A *= Sq(D(1) + (0.37464 + 1.54226*w - 0.26992*w*w)*(D(1) - sqrt(Tr))); }
				   break;
	}


	switch (FModel)
	{
	case FI_SC:    Fe = Tau;   Fm = 0.5;   break;
	case FI_KUPER: Fe = 1.0;   Fm = 0.5;     break;
	}

	string FileName = Name(MModel);
	FileName = "../data/CoCurve_" + FileName.substr(3, 3) + "_2.txt";

	ifstream File;
	File.open(FileName);
	if (!File.is_open())
	{
		cout << "open file error:  " << FileName << endl;
		return;
	}

	int TemNum = 0;
	char   Buffer[512];
	double TemDen[1000][3];
	istringstream Iss;

	while (!File.eof())
	{
		double T0, DenGas, DenLiquid;
		File.getline(Buffer, 512);
		Iss.clear();  Iss.str(Buffer);
		Iss >> T0 >> DenGas >> DenLiquid;

		if (T0 > 0 && T0 < 1.1)
		{
			TemDen[TemNum][0] = T0;
			TemDen[TemNum][1] = DenGas;
			TemDen[TemNum][2] = DenLiquid;
			++TemNum;
		}
	}
	File.close();


	TemDen[TemNum][0] = 1.12;

	for (int n = 0; n < TemNum; ++n)
	{
		if (Eq(Tr, TemDen[n][0]))
		{
			T = TemDen[n][0] * Tc;
			DenG = TemDen[n][1] * Rc;
			DenL = TemDen[n][2] * Rc;
			break;
		}
	}

}

//*************************************************************************************************
void LocalCollideMrt(const GridData & t, double * Df)
{

	static const double M[DQ][DQ] = {
		1, 1, 1, 1, 1, 1, 1, 1, 1,
		-4,-1,-1,-1,-1, 2, 2, 2, 2,
		4,-2,-2,-2,-2, 1, 1, 1, 1,
		0, 1, 0,-1, 0, 1,-1,-1, 1,
		0,-2, 0, 2, 0, 1,-1,-1, 1,
		0, 0, 1, 0,-1, 1, 1,-1,-1,
		0, 0,-2, 0, 2, 1, 1,-1,-1,
		0, 1,-1, 1,-1, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 1,-1, 1,-1
	};

	static const double R[DQ][DQ] = {
		1. / 9, -1. / 9,   1. / 9,      0,      0,     0,     0,     0,     0,
		1. / 9, -1. / 36, -1. / 18,  1. / 6,  -1. / 6,     0,     0,  1. / 4,     0,
		1. / 9, -1. / 36, -1. / 18,     0,      0,  1. / 6, -1. / 6, -1. / 4,     0,
		1. / 9, -1. / 36, -1. / 18, -1. / 6,   1. / 6,     0,     0,  1. / 4,     0,
		1. / 9, -1. / 36, -1. / 18,     0,      0, -1. / 6,  1. / 6, -1. / 4,     0,
		1. / 9,  1. / 18,  1. / 36,  1. / 6,  1. / 12,  1. / 6,  1. / 12,    0,  1. / 4,
		1. / 9,  1. / 18,  1. / 36, -1. / 6, -1. / 12,  1. / 6,  1. / 12,    0, -1. / 4,
		1. / 9,  1. / 18,  1. / 36, -1. / 6, -1. / 12, -1. / 6, -1. / 12,    0,  1. / 4,
		1. / 9,  1. / 18,  1. / 36,  1. / 6,  1. / 12, -1. / 6, -1. / 12,    0, -1. / 4
	};

	const double Ts[DQ] = { 0, 1.64, 1.54, 0, 1.54, 0, 1.54, ReTau, ReTau };

	static double Meq[DQ], Mft[DQ], Mf[DQ];

	Meq[0] = t.Den;
	Meq[1] = t.Den * (-2.0 + 3.0*(Sq(t.Vx) + Sq(t.Vy)));
	Meq[2] = t.Den * (1.00 - 3.0*(Sq(t.Vx) + Sq(t.Vy)));
	Meq[3] = t.Den * t.Vx;
	Meq[4] = t.Den * t.Vx * -1;
	Meq[5] = t.Den * t.Vy;
	Meq[6] = t.Den * t.Vy * -1;
	Meq[7] = t.Den * (Sq(t.Vx) - Sq(t.Vy));
	Meq[8] = t.Den * t.Vx * t.Vy; //*/

	for (int f = 0; f < DQ; ++f)
	{
		Mf[f] = M[f][0] * Df[0] + M[f][1] * Df[1] + M[f][2] * Df[2] + M[f][3] * Df[3] + M[f][4] * Df[4] + M[f][5] * Df[5] + M[f][6] * Df[6] + M[f][7] * Df[7] + M[f][8] * Df[8];
	}

	for (int f = 0; f < DQ; ++f)
	{
		Mf[f] = Mf[f] - Ts[f] * (Mf[f] - Meq[f]);
	}

	for (int f = 0; f < DQ; ++f)
	{
		Df[f] = R[f][0] * Mf[0] + R[f][1] * Mf[1] + R[f][2] * Mf[2] + R[f][3] * Mf[3] + R[f][4] * Mf[4] + R[f][5] * Mf[5] + R[f][6] * Mf[6] + R[f][7] * Mf[7] + R[f][8] * Mf[8];
	}
}
//*************************************************************************************************
void GlobalCollide()
{
#pragma omp parallel for 
	FOR_iDX_jDY_Fluid
	{
		GridData & t = Grid[i][j];

	double   * Df = Dist[i][j];
		t.EVx = t.Vx + t.Fx / t.Den * Fe;
		t.EVy = t.Vy + t.Fy / t.Den * Fe;
		t.MVx = t.Vx + t.Fx / t.Den * Fm;
		t.MVy = t.Vy + t.Fy / t.Den * Fm;

		LocalCollideMrt(t, Df);
		for (int f = 0; f < DQ; ++f)
		{
			int ii = i + Ex[f];	if (ii < 0) ii += DX;	else if (ii >= DX) ii -= DX;
			int jj = j + Ey[f];	if (jj < 0) jj += DY;	else if (jj >= DY) jj -= DY;
			switch (FModel)
			{
			case FI_SC: {
				(Grid[ii][jj].Type == FLUID ? Temp[ii][jj][f] : Temp[i][j][Re[f]]) = Df[f];
			}; break;
			case FI_KUPER: {

				Df[f] += Feq(f, t.Den, t.EVx, t.EVy) - Feq(f, t.Den, t.Vx, t.Vy);
				(Grid[ii][jj].Type == FLUID ? Temp[ii][jj][f] : Temp[i][j][Re[f]]) = Df[f];
			}break;
			default:
				cout << "unknown forcing item:  " << FModel << endl;
			}
		}
		}

	double(*p)[DY][DQ];
	p = Dist;   Dist = Temp;   Temp = p;

	}
//*************************************************************************************************
void MacroCalculate(void)
{
	double TmpMass = 0;
#pragma omp parallel for reduction(+:TmpMass) 
	FOR_iDX_jDY_Fluid
	{
		GridData & t = Grid[i][j];
		double  *  Df = Dist[i][j];

		t.Den = Df[0] + Df[1] + Df[2] + Df[3] + Df[4] + Df[5] + Df[6] + Df[7] + Df[8];
		t.Vx = (Df[1] + Df[5] + Df[8] - Df[3] - Df[6] - Df[7]) / t.Den;
		t.Vy = (Df[2] + Df[5] + Df[6] - Df[4] - Df[7] - Df[8]) / t.Den;

		TmpMass += t.Den;

		if (t.Den != t.Den || t.Den < 0 || fabs(t.Den - Density)>10)
		{
			cout << "Density: " << NowStep << ",  (" << i << "," << j << ")   " << t.Den << endl;
		}
	}
	Mass = TmpMass;

	//===============================================================

	if (NowStep%min(ShowStep, SaveStep) == 0)
	{
		MaxSpeed = 0;
		FOR_iDX_jDY_Fluid
		{
			double Mod = Grid[i][j].MVy;
			if (Mod > MaxSpeed)   MaxSpeed = Mod;
		}

	}
}

//*************************************************************************************************
void NonidealForce()
{

#pragma omp parallel for
	FOR_iDX_jDY_Fluid
	{
	GridData & t = Grid[i][j];
	Define_ij6;
	t.Dx = DerivativeX6(Den);
	t.Dy = DerivativeY6(Den);
	t.Dd = DerivativeD6(Den);

		switch (MModel)
		{
			case MP_P0VDW:
				t.P0 = T * t.Den / (D(1) - B * t.Den) - A * Sq(t.Den);
				break;

			case MP_P0CSE:
				t.P0 = T * t.Den*((D(1) + (B*t.Den) + Sq(B*t.Den) - Cu(B*t.Den)) / Cu(D(1) - B * t.Den)) - A * Sq(t.Den);
				break;

			case MP_P0PRM:
			case MP_P0PRW:
				t.P0 = T * t.Den / (D(1) - B * t.Den) - A * Sq(t.Den) / (D(1) + B * t.Den * 2 - Sq(B*t.Den));
				break;

			case MP_P0RKE:
			case MP_P0RKS:
				t.P0 = T * t.Den / (D(1) - B * t.Den) - A * Sq(t.Den) / (D(1) + B * t.Den);
				break;

			default:
				cout << "unknown multiphase model: " << MModel << endl;
			}

			t.P0 = Sq(Kr) * t.P0;
	}

#pragma omp parallel for
		FOR_iDX_jDY_Fluid
	{
	GridData & t = Grid[i][j];
		Define_ij6;
	double Ddx = DerivativeX6(Dd);
	double Ddy = DerivativeY6(Dd);

	t.Fx = -DerivativeX6(P0) + K * t.Den*Ddx + t.Dx / 3;
	t.Fy = -DerivativeY6(P0) + K * t.Den*Ddy + t.Dy / 3;
	}
}

//*************************************************************************************************
void ShowData()
{
	cout << setw(10) << NowStep << "    " << setiosflags(ios::fixed) << setprecision(12) << Mass << setprecision(5) << setw(12) << Grid[DX / 2][0].Den << setw(12) << Grid[DX / 2][DY / 2].Den << setw(12) << MaxSpeed << setw(10) << StepTime << setw(10) << Tr << endl;
	
}

void SaveData()
{
	char FileTecplot[512];
	sprintf(FileTecplot, "../results/Tecplot_%3.1f_%d.dat", Tau, NowStep);  
	File.open(FileTecplot);
	File << "TITLE = SPLASHING" << endl;
	File << "VARIABLES =  X,   Y,   Density,   U,   V" << endl;
	File << "ZONE  I=" << DX << "  J=" << DY - 2 << "  F=POINT" << endl;

	for (int j = 1; j < DY - 1; ++j)
	{
		for (int i = 0; i < DX; ++i)
		{
			GridData & t = Grid[(i) % DX][j];
			File << D(i) << "  " << j - 1 << "  " << setprecision(6) << (t.Type == FLUID ? t.Den : -1) << "  " << t.MVx << "  " << t.MVy << endl;
		}
	}
}


