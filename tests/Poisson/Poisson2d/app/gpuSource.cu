template <typename T>  __global__  void kernelgpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		f[0*ng+i] = sin(xdg1*3.141592653589793)*sin(xdg2*3.141592653589793)*1.973920880217872E1;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuSource<<<gridDim, blockDim>>>(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
