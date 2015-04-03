
typedef struct {
  int nlon;
  int nlat;
  int ntr;
  int poldim;
  int ifax[10];
  double *trig;
  double *poli;
  double *pold;
  double *pol2;    /* only for uv2dv  */
  double *pol3;    /* only for uv2dv  */
  double *coslat;  /* only for scaluv with uv2dv */
  double *rcoslat; /* only for scaluv with dv2uv */
}
SPTRANS;

typedef struct {
  int ntr;
  int fdim;
  double *f1;
  double *f2;
}
DVTRANS;

void dv2ps(double *div, double *pot, int nlev, int ntr);

SPTRANS *sptrans_new(int nlon, int nlat, int ntr, int flag);
void sptrans_delete(SPTRANS *sptrans);

DVTRANS *dvtrans_new(int ntr);
void dvtrans_delete(DVTRANS *dvtrans);

void trans_uv2dv(SPTRANS *sptrans, int nlev,
		 int gridID1, double *gu, double *gv,
		 int gridID2, double *sd, double *svo);

void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev,
		 int gridID1, double *sd, double *svo,
		 int gridID2, double *gu, double *gv);

void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);

void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void speccut(int gridIDin, double *arrayIn, double *arrayOut, int *waves);

void sp2fc(double *sa, double *fa, double *poli, int nlev, int nlat, int nfc, int nt);
void fc2sp(double *fa, double *sa, double *poli, int nlev, int nlat, int nfc, int nt);

void fc2gp(double *trig, int *ifax, double *fc, double *gp, int nlat, int nlon, int nlev, int nfc);
void gp2fc(double *trig, int *ifax, double *gp, double *fc, int nlat, int nlon, int nlev, int nfc);

void sp2sp(double *arrayIn, int ntrIn, double *arrayOut, int ntrOut);
void spcut(double *arrayIn, double *arrayOut, int ntr, int waves[]);

void phcs(double *pnm, double *hnm, int waves, double pmu,
	  double *ztemp1, double *ztemp2);
void jspleg1(double *pleg, double plat, int ntr, double *work);

void fft_set(double *trigs, int *ifax, int n);
