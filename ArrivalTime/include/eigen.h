#ifdef __cpluscplus
extern "C" {
#endif
int bakvec ( int n, float t[], float e[], int m, float z[] );
void balbak ( int n, int low, int igh, float scale[], int m, float z[] );
void bandr ( int n, int mb, float a[], float d[], float e[], float e2[], 
  int matz, float z[] );
void cbabk2 ( int n, int low, int igh, float scale[], int m, float zr[], 
  float zi[] );
void csroot ( float xr, float xi, float *yr, float *yi );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
float pythagoras ( float a, float b );
float r8_abs ( float x );
float r8_epsilon ( void );
float r8_max ( float x, float y );
float r8_min ( float x, float y );
float r8_sign ( float x );
void r8mat_identity  ( int n, float a[] );
float *r8mat_mm_new ( int n1, int n2, int n3, float a[], float b[] );
void r8mat_print ( int m, int n, float a[], char *title );
void r8mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
float *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8vec_print ( int n, float a[], char *title );
int rs ( int n, float a[], float w[], int matz, float z[] );
int rsb ( int n, int mb, float a[], float w[], int matz, float z[] );
void timestamp ( void );
int tql22 ( int n, float d[], float e[], float z[] );
int tqlrat ( int n, float w[], float fv2[] );
void tred1 ( int n, float a[], float w[], float fv1[], float fv2[] );
void tred2 ( int n, float a[], float d[], float e[], float z[] );
#ifdef _cplusplus
}
#endif

