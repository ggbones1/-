#pragma once
#pragma warning( disable : 4018 )  
#pragma warning( disable : 4267 )  
#pragma warning( disable : 4786 )  
#pragma warning( disable : 4996 )  

//*************************************************************************************************
const int DX = 401;
const int DY = 401;

#define LBEMRT
//*************************************************************************************************
const char FLUID = 1;
const char BOUNDARY = 2;
const char UNKNOWN = 97;

const int MP_P0VDW = 116;
const int MP_P0PRW = 117;
const int MP_P0RKE = 118;
const int MP_P0RKS = 119;
const int MP_P0CSE = 120;
const int MP_P0PRM = 121;

const int FI_SC = 200;
const int FI_KUPER = 201;

inline const char * Name(const int Model)
{
	switch (Model)
	{
	case MP_P0VDW:	return "P0_VDW";
	case MP_P0PRW:	return "P0_PRW";
	case MP_P0RKE:  return "P0_RKE";
	case MP_P0RKS:  return "P0_RKS";
	case MP_P0CSE:  return "P0_CSE";
	case MP_P0PRM:	return "P0_PRM";
	case FI_SC:     return "SC";
	case FI_KUPER:  return "KS";
	}
	return "Unknown";
}

//*************************************************************************************************
const int DQ = 9;
const int Ex[DQ] = { 0,  1,  0, -1,  0,  1, -1, -1,  1 };
const int Ey[DQ] = { 0,  0,  1,  0, -1,  1,  1, -1, -1 };
const int Re[DQ] = { 0,  3,  4,  1,  2,  7,  8,  5,  6 };
const double Gamma[DQ] = { 0, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 12, 1.0 / 12, 1.0 / 12, 1.0 / 12 };
const double Alpha[DQ] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };
const double PI = 3.14159265358979323846264338327950288;
const double S2 = 1.4142135623730950488016887242097;

//*************************************************************************************************
#define D(x)  (double(x))
#define FOR_iDX_jDY	 for(int i=0; i<DX;++i) for(int j=0;j<DY;++j)
#define FOR_iDX_jDY_Fluid  for(int i=0; i<DX;++i) for(int j=0;j<DY;++j) if(Grid[i][j].Type == FLUID)

#define Define_i1i2j1j2   int i1 = (i<DX-1?i+1:0), i2 = (i>0?i-1:DX-1),  j1 = (j<DY-1?j+1:0), j2 = (j>0?j-1:DY-1)
#define DerivativeX( x )  (( Grid[i1][j].x-Grid[i2][j].x)/3 + (Grid[i1][j1].x-Grid[i2][j2].x)/12 + (Grid[i1][j2].x-Grid[i2][j1].x)/12)
#define DerivativeY( x )  (( Grid[i][j1].x-Grid[i][j2].x)/3 + (Grid[i1][j1].x-Grid[i2][j2].x)/12 + (Grid[i2][j1].x-Grid[i1][j2].x)/12)
#define DerivativeD( x )  (( Grid[i][j1].x + Grid[i][j2].x  +  Grid[i1][j].x +Grid[i2][j].x)*2/3 + (Grid[i1][j1].x+Grid[i1][j2].x + Grid[i2][j1].x + Grid[i2][j2].x)/3 - Grid[i][j].x*4)

#define Define_ij5  int i2 = (i >0?i-1 :DX-1), i4 = (i <DX-1?i+1 :0), j2 = (j >0?j-1 :DY-1), j4 = (j <DY-1?j+1 :0); \
				    int i1 = (i2>0?i2-1:DX-1), i5 = (i4<DX-1?i4+1:0), j1 = (j2>0?j2-1:DY-1), j5 = (j4<DY-1?j4+1:0) 

#define DerivativeX5( x )  (  \
(Grid[i1][j ].x - Grid[i2][j ].x * 8 + Grid[i4][j ].x * 8 - Grid[i5][j ].x) / 18 +  \
(Grid[i1][j1].x - Grid[i2][j2].x * 8 + Grid[i4][j4].x * 8 - Grid[i5][j5].x) / 72 +  \
(Grid[i1][j5].x - Grid[i2][j4].x * 8 + Grid[i4][j2].x * 8 - Grid[i5][j1].x) / 72 )

#define DerivativeY5( x )  (  \
(Grid[i ][j1].x - Grid[i ][j2].x * 8 + Grid[i ][j4].x * 8 - Grid[i ][j5].x) / 18 +  \
(Grid[i1][j1].x - Grid[i2][j2].x * 8 + Grid[i4][j4].x * 8 - Grid[i5][j5].x) / 72 +  \
(Grid[i5][j1].x - Grid[i4][j2].x * 8 + Grid[i2][j4].x * 8 - Grid[i1][j5].x) / 72 )

#define DerivativeD5( x )  (  \
(-Grid[i1][j].x + Grid[i2][j].x * 16 - Grid[i][j].x * 30 + Grid[i4][j].x * 16 - Grid[i5][j].x) / 12 + \
(-Grid[i][j1].x + Grid[i][j2].x * 16 - Grid[i][j].x * 30 + Grid[i][j4].x * 16 - Grid[i][j5].x) / 12 )
//*************************************************************************************************
#define Define_ij6 int i3=(i>0?i-1:DX-1),i5=(i<DX-1?i+1:0),j3=(j>0?j-1:DY-1),j5=(j<DY-1?j+1:0);\
				   int i2=(i3>0?i3-1:DX-1),i6=(i5<DX-1?i5+1:0),j2=(j3>0?j3-1:DY-1),j6=(j5<DY-1?j5+1:0);\
				   int i1=(i2>0?i2-1:DX-1),i7=(i6<DX-1?i6+1:0),j1=(j2>0?j2-1:DY-1),j7=(j6<DY-1?j6+1:0)
#define DerivativeX6(x)  (  \
		(Grid[i7][j ].x - Grid[i1][j ].x -9*Grid[i6][j].x+ 9*Grid[i2][j].x + 45*Grid[i5][j ].x-45* Grid[i3][j ].x)/ 90 +  \
		(Grid[i7][j7].x - Grid[i1][j1].x -9*Grid[i6][j6].x+ 9*Grid[i2][j2].x+ 45*Grid[i5][j5].x-45*Grid[i3][j3].x) / 360 +  \
		(Grid[i7][j1].x - Grid[i1][j7].x -9*Grid[i6][j2].x +9* Grid[i2][j6].x  +45* Grid[i5][j3].x-45*Grid[i3][j5].x) / 360)
#define DerivativeY6(x)  (  \
		(Grid[i][j7].x - Grid[i][j1].x -9*Grid[i][j6].x+ 9*Grid[i][j2].x + 45*Grid[i][j5 ].x-45* Grid[i][j3 ].x)/ 90 +  \
		(Grid[i7][j7].x - Grid[i1][j1].x -9*Grid[i6][j6].x+ 9*Grid[i2][j2].x+ 45*Grid[i5][j5].x-45*Grid[i3][j3].x) / 360 +  \
		(Grid[i1][j7].x - Grid[i7][j1].x -9*Grid[i2][j6].x +9* Grid[i6][j2].x  +45* Grid[i3][j5].x-45*Grid[i5][j3].x) / 360)
#define DerivativeD6( x )  (  (\
		(2*Grid[i1][j].x  -27*Grid[i2][j].x +270*Grid[i3][j].x+270*Grid[i5][j].x- 27*Grid[i6][j].x+2*Grid[i7][j].x ) *4 + \
		(2*Grid[i1][j1].x  -27*Grid[i2][j2].x +270*Grid[i3][j3].x+270*Grid[i5][j5].x- 27*Grid[i6][j6].x+2*Grid[i7][j7].x)  +\
		(2*Grid[i1][j7].x  -27*Grid[i2][j6].x +270*Grid[i3][j5].x+270*Grid[i5][j3].x- 27*Grid[i6][j2].x+2*Grid[i7][j1].x )  +\
		( 2*Grid[i][j1].x  - 27*Grid[i][j2].x +270*Grid[i][j3].x +270*Grid[i][j5].x-27*Grid[i][j6].x+2*Grid[i][j7].x  ) *4+\
		(2 * Grid[i1][j1].x - 27 * Grid[i2][j2].x + 270 * Grid[i3][j3].x + 270 * Grid[i5][j5].x - 27 * Grid[i6][j6].x + 2 * Grid[i7][j7].x ) +\
		( 2*Grid[i7][j1].x  - 27*Grid[i6][j2].x +270*Grid[i5][j3].x +270*Grid[i3][j5].x-27*Grid[i2][j6].x+2*Grid[i1][j7].x  )-5880*Grid[i][j].x )/1080)

#define Define_ij4  int i2 = (i >0?i-1 :DX-1), i4 = (i <DX-1?i+1 :0), j2 = (j >0?j-1 :DY-1), j4 = (j <DY-1?j+1 :0); \
				    int i1 = (i2>0?i2-1:DX-1), i5 = (i4<DX-1?i4+1:0), j1 = (j2>0?j2-1:DY-1), j5 = (j4<DY-1?j4+1:0)
#define DerivativeX4( x ) (  \
		(Grid[i1][j ].x - Grid[i2][j ].x * 8 + Grid[i4][j ].x * 8 - Grid[i5][j ].x)/ 18 +  \
		(Grid[i1][j1].x - Grid[i2][j2].x * 8 + Grid[i4][j4].x * 8 - Grid[i5][j5].x) / 72 +  \
		(Grid[i1][j5].x - Grid[i2][j4].x * 8 + Grid[i4][j2].x * 8 - Grid[i5][j1].x) / 72 )
#define DerivativeY4( x ) (  \
		(Grid[i ][j1].x - Grid[i ][j2].x * 8 + Grid[i ][j4].x * 8 - Grid[i ][j5].x) /18 +  \
		(Grid[i1][j1].x - Grid[i2][j2].x * 8 + Grid[i4][j4].x * 8 - Grid[i5][j5].x) / 72 +  \
		(Grid[i5][j1].x - Grid[i4][j2].x * 8 + Grid[i2][j4].x * 8 - Grid[i1][j5].x) / 72 )
#define DerivativeD4( x )  (( \
(-Grid[i1][j ].x + Grid[i2][j ].x * 16  + Grid[i4][j ].x * 16 - Grid[i5][j ].x) * 4 +  \
(-Grid[i1][j5].x + Grid[i2][j4].x * 16  + Grid[i4][j2].x * 16 - Grid[i5][j1].x)     +  \
(-Grid[i1][j1].x + Grid[i2][j2].x * 16  + Grid[i4][j4].x * 16 - Grid[i5][j5].x)     +  \
(-Grid[i ][j1].x + Grid[i][j2 ].x * 16  + Grid[i][j4 ].x * 16 - Grid[i ][j5].x) * 4 +  \
(-Grid[i5][j1].x + Grid[i4][j2].x * 16  + Grid[i2][j4].x * 16 - Grid[i1][j5].x)     +  \
(-Grid[i1][j1].x + Grid[i2][j2].x * 16  + Grid[i4][j4].x * 16 - Grid[i5][j5].x) - Grid[i][j].x * 360) / 72)

#define Define_ij2   int i1 = (i<DX-1?i+1:0), i2 = (i>0?i-1:DX-1),  j1 = (j<DY-1?j+1:0), j2 = (j>0?j-1:DY-1)
#define DerivativeX2( x )  (( Grid[i1][j].x-Grid[i2][j].x)/3 + (Grid[i1][j1].x-Grid[i2][j2].x)/12 + (Grid[i1][j2].x-Grid[i2][j1].x)/12)
#define DerivativeY2( x )  (( Grid[i][j1].x-Grid[i][j2].x)/3 + (Grid[i1][j1].x-Grid[i2][j2].x)/12 + (Grid[i2][j1].x-Grid[i1][j2].x)/12)
#define DerivativeD2( x )  (( Grid[i][j1].x+Grid[i][j2].x    +  Grid[i1][j].x +Grid[i2][j].x)*2/3 + (Grid[i1][j1].x+Grid[i1][j2].x + Grid[i2][j1].x + Grid[i2][j2].x)/3 - Grid[i][j].x*4)

template<class T> inline T Sq(T x)
{
	return x * x;
}

template<class T> inline T Cu(T x)
{
	return x * x*x;
}


inline bool Eq(double x, double y, double d = 1E-12)
{
	return (x > y ? x - y : y - x) < d;
}

inline double Feq(int f, double Den, double Vx, double Vy)
{
	double DotMet = Vx * Ex[f] + Vy * Ey[f];
	return Den * Alpha[f] * (1.0 + 3.0*DotMet + 4.5*DotMet*DotMet - 1.5*(Vx*Vx + Vy * Vy));
}

inline int iModx(int iData)
{
	return (iData + DX) % DX;
}
inline int iMody(int iData)
{
	return (iData + DY) % DY;
}
