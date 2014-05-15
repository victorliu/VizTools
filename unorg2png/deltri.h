typedef struct delaunay_triangulation *deltri;

deltri deltri_new(int n, double *v);
int DT_interpolate(deltri T, const double xy[2], int *h, double *val); // returns 0 if within the mesh
void DT_destroy(deltri T);
