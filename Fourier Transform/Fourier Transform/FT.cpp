#include "FT.h"

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i<M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j<N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	//initialization
	int M = h, N = w;
	double * xreal;double *ximg;
	double ** resultReal;
	double ** resultImg;
	xreal = new double[M];
	ximg = new double[M];
	resultReal = new double*[M];
	resultImg = new double*[M];
	for (int index = 0; index < M; index++) {
		resultReal[index] = new double[N];
		resultImg[index] = new double[N];
	}
	for (int indexa = 0; indexa < M; indexa++) {
		for (int indexb = 0; indexb < M; indexb++)
		{
			resultReal[indexa][indexb] = InputImage[indexb][indexa];
			resultImg[indexa][indexb] = 0;
		}
	}
	//do 1d fft to every column
	for (int indexa = 0; indexa < M; indexa++) {
		for (int indexb = 0; indexb < M; indexb++) {
			xreal[indexb] = resultReal[indexb][indexa];
			ximg[indexb] = resultImg[indexb][indexa];
		}
		FFT(resultReal,resultImg,InputImage,xreal,ximg,M,N,indexa,indexa);
		for (int indexb = 0; indexb < M; indexb++) {
			resultReal[indexb][indexa] = xreal[indexb];
			resultImg[indexb][indexa] = ximg[indexb];
		}
	}
	//do 1d fft to every row
	for (int indexa = 0; indexa < M; indexa++) {
		for (int indexb = 0; indexb < M; indexb++) {
			xreal[indexb] = resultReal[indexa][indexb];
			ximg[indexb] = resultImg[indexa][indexb];
		}
		FFT(resultReal, resultImg, InputImage, xreal, ximg, M, N, indexa, indexa);
		for (int indexb = 0; indexb < M; indexb++) {
			resultReal[indexa][indexb] = xreal[indexb];
			resultImg[indexa][indexb] = ximg[indexb];
		}
	}
	//send transformed information to output
	for (int indexa = 0; indexa < M; indexa++) {
		for (int indexb = 0; indexb < M; indexb++) {
			FreqImag[indexa][indexb] = resultImg[indexa][indexb];
			FreqReal[indexa][indexb] = resultReal[indexa][indexb];
			OutputImage[indexa][indexb] = sqrt(resultReal[indexa][indexb] * resultReal[indexa][indexb] + resultImg[indexa][indexb] * resultImg[indexa][indexb]) * 255;
		}
	}

	for (int index = 0; index < M; index++) {
		delete[] resultReal[index];
		delete[] resultImg[index];
	}
	delete[] resultReal;
	delete[] resultImg;

	delete[] xreal;
	delete[] ximg;
}

void FT::FFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage,double * xreal,double *ximg, int h, int w, int u, int v)
{
	int M = h;
	int N = w;

	int SIZE = M;
	for (int indexa = 1, indexb = 0; indexa < SIZE; indexa++) {
		for (int k = SIZE >> 1; !((indexb ^= k)&k); k >>= 1);

		if (indexa > indexb) {
			double tempr = xreal[indexa], tempi = ximg[indexa];
			xreal[indexa] = xreal[indexb]; ximg[indexa] = ximg[indexb];
			xreal[indexb] = tempr; ximg[indexb] = tempi;
		}
	}

	for (int k = 2; k <= SIZE; k <<= 1) {

		double w = -2 * 3.14159 / (1.0*k);
		double dsitar = cos(w), dsitai = sin(w);

		for (int t = 0; t < SIZE; t += k) {

			double sitar = 1, sitai = 0;

			for (int tt = t; tt < t + k / 2; tt++) {
				double ar = xreal[tt];
				double ai = ximg[tt];
				double br = xreal[tt + k / 2] * sitar - ximg[tt + k / 2] * sitai;
				double bi = xreal[tt + k / 2] * sitai + ximg[tt + k / 2] * sitar;
				xreal[tt] = ar + br;
				ximg[tt] = ai + bi;
				xreal[tt + k / 2] = ar - br;
				ximg[tt + k / 2] = ai - bi;
				double tempr = sitar, tempi = sitai;
				sitar = tempr*dsitar - tempi*dsitai;
				sitai = tempr*dsitai + tempi*dsitar;
			}
		}
	}
	for (int indexa = 0; indexa < SIZE; indexa++) {
		xreal[indexa] /= SIZE;
		ximg[indexa] /= SIZE;
	}

}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	std::cout << "hello,world" << std::endl;



	int SIZE = h;
	double ** resultReal;
	double ** resultImag;
	resultReal = new double*[SIZE];
	resultImag = new double*[SIZE];
	for (int index = 0; index < SIZE; index++) {
		resultReal[index] = new double[SIZE];
		resultImag[index] = new double[SIZE];
	}

	double * xreal = new double[SIZE];
	double * ximag = new double[SIZE];

	for (int indexa = 0; indexa < SIZE; indexa++) {
		for (int indexb = 0; indexb < SIZE; indexb++) {
			resultReal[indexa][indexb] = FreqReal[indexb][indexa];
			resultImag[indexa][indexb] = FreqImag[indexb][indexa];
		}
	}



	for (int indexa = 0; indexa < SIZE; indexa++) {
		for (int indexb = 0; indexb < SIZE; indexb++) {
			xreal[indexb] = resultReal[indexa][indexb] * SIZE;
			ximag[indexb] = resultImag[indexa][indexb] * SIZE;
		}
		InverseFFT(FreqReal,FreqImag,FreqReal,FreqImag,xreal,ximag,h,w,indexa,indexa);
		for (int indexb = 0; indexb < SIZE; indexb++) {
			resultReal[indexa][indexb] = xreal[indexb];
			resultImag[indexa][indexb] = ximag[indexb];
		}
	}

	for (int indexa = 0; indexa < SIZE; indexa++) {
		for (int indexb = 0; indexb < SIZE; indexb++) {
			xreal[indexb] = resultReal[indexb][indexa]*SIZE;
			ximag[indexb] = resultImag[indexb][indexa]*SIZE;
		}

		InverseFFT(FreqReal, FreqImag, FreqReal, FreqImag, xreal, ximag, h, w, indexa, indexa);
		for (int indexb = 0; indexb < SIZE; indexb++) {
			resultReal[indexb][indexa] = xreal[indexb];
			resultImag[indexb][indexa] = ximag[indexb];
		}
	}



	double flag = -100;
	for (int indexa = 0; indexa < SIZE; indexa++) {
		for (int indexb = 0; indexb < SIZE; indexb++) {
			OutputImage[indexa][indexb] = sqrt(resultReal[indexa][indexb] * resultReal[indexa][indexb] + resultImag[indexa][indexb] * resultImag[indexa][indexb]);
			if (flag < OutputImage[indexa][indexb])	flag = OutputImage[indexa][indexb];
		}
	}
	/*
	for (int index = 0; index < w; index++) {
		for (int indexa = 0; indexa < w; indexa++) {
			std::cout << OutputImage[index][indexa]<<" ";
		}
		std::cout << std::endl;
	}*/
	std::cout << flag << std::endl;
	flag /= 255;
	flag = 1 / flag;

	for (int indexa = 0; indexa < SIZE; indexa++) {
		for (int indexb = 0; indexb < SIZE; indexb++) {
			FreqImag[indexa][indexb] = resultImag[indexa][indexb];
			FreqReal[indexa][indexb] = resultReal[indexa][indexb];
			OutputImage[indexa][indexb] *= flag;
		}
	}
	//TO-DO
	for (int i = SIZE - 1; i >= SIZE / 2; i--)
	{
		for (int j = SIZE - 1; j >= 0; j--)
		{
			double temp = OutputImage[i][j];
			OutputImage[i][j] = OutputImage[SIZE - i][j];
			OutputImage[SIZE - i][j] = temp;
		}
	}
	for (int i = SIZE - 1; i >= 0; i--)
	{
		for (int j = SIZE - 1; j >= SIZE / 2; j--)
		{
			double temp = OutputImage[i][j];
			OutputImage[i][j] = OutputImage[i][SIZE - j];
			OutputImage[i][SIZE - j] = temp;
		}
	}


	for (int index = 0; index < SIZE; index++) {
		delete[] resultReal[index];
		delete[] resultImag[index];
	}
	delete[] resultReal;
	delete[] resultImag;
	delete[] xreal;
	delete[] ximag;
}

void FT::InverseFFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag,double * xreal,double * ximag, int h, int w, int x, int y)
{
	int SIZE = w;
	for (int indexa = 1, indexb = 0; indexa < SIZE; indexa++) {
		for (int k = SIZE >> 1; !((indexb ^= k)&k); k >>= 1);

		if (indexa > indexb) {
			double tempr = xreal[indexa], tempi = ximag[indexa];
			xreal[indexa] = xreal[indexb]; ximag[indexa] = ximag[indexb];
			xreal[indexb] = tempr; ximag[indexb] = tempi;
		}
	}

	for (int k = 2; k <= SIZE; k <<= 1) {
		double w = -2 * 3.14159 / (1.0*k);
		double dsitar = cos(w), dsitai = sin(w);

		for (int t = 0; t < SIZE; t += k) {

			double sitar = 1, sitai = 0;

			for (int tt = t; tt < t + k / 2; tt++) {
				double ar = xreal[tt];
				double ai = ximag[tt];
				double br = xreal[tt + k / 2] * sitar - ximag[tt + k / 2] * sitai;
				double bi = xreal[tt + k / 2] * sitai + ximag[tt + k / 2] * sitar;
				xreal[tt] = ar + br;
				ximag[tt] = ai + bi;
				xreal[tt + k / 2] = ar - br;
				ximag[tt + k / 2] = ai - bi;
				double tempr = sitar, tempi = sitai;
				sitar = tempr*dsitar - tempi*dsitai;
				sitai = tempr*dsitai + tempi*dsitar;
			}
		}
	}
}


//team 11's code
void swap(std::complex<double> *a, std::complex<double> *b)
{
	std::complex<double> temp = *a;
	*a = *b;
	*b = temp;
}


void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w, int dir)
{
	int M = h;
	int N = w;
	printf("M = %d,N = %d\n", M, N);

	std::complex<double>** out = new std::complex<double>*[M];
	for (int newcnt = 0; newcnt<M; newcnt++)
	{
		out[newcnt] = new std::complex<double>[N]; // 傅立葉頻率陣列
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double a = FreqReal[j][i], b = FreqImag[j][i];
			if (dir == 1)
			{
				a = InputImage[j][i];
				b = 0;
			}
			std::complex<double> temp(a, b);
			out[i][j] = temp;

		}
	}
	if (dir == 1)
		for (int j = 0; j <N; ++j)
		{
			std::complex<double> *input = new std::complex<double>[M];
			for (int i = 0; i < M; ++i)
				input[i] = out[i][j];
			FFT(FreqReal[0], FreqImag[0], input, M, M, dir);
			for (int i = 0; i < M; ++i)
				out[i][j] = input[i];
			free(input);
		}
	for (int i = 0; i < M; ++i)
	{
		std::complex<double> *input = new std::complex<double>[N];
		for (int j = 0; j < N; ++j)
		{
			if (dir == 1)
				input[j] = out[i][j];
			else
			{
				std::complex<double> temp(out[i][j].real()*N, out[i][j].imag()*N);
				input[j] = temp;
			}
		}
		FFT(FreqReal[0], FreqImag[0], input, N, N, dir);
		for (int j = 0; j < N; ++j)
			out[i][j] = input[j];
		free(input);
	}
	if (dir == -1)
		for (int j = 0; j <N; ++j)
		{
			std::complex<double> *input = new std::complex<double>[M];
			for (int i = 0; i < M; ++i)
			{
				std::complex<double> temp(out[i][j].real()*N, out[i][j].imag()*N);
				input[i] = temp;
			}
			FFT(FreqReal[0], FreqImag[0], input, M, M, dir);
			for (int i = 0; i < M; ++i)
				out[i][j] = input[i];
			free(input);
		}
	double max_s = -100;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			FreqImag[i][j] = out[i][j].imag();
			FreqReal[i][j] = out[i][j].real();
			if (dir == 1)
				OutputImage[i][j] = sqrt(pow(out[i][j].real(), 2.0) + pow(out[i][j].imag(), 2.0)) * 255;
			else
			{
				OutputImage[i][j] = sqrt(pow(out[i][j].real(), 2.0) + pow(out[i][j].imag(), 2.0));
			}
			max_s < OutputImage[i][j] ? max_s = OutputImage[i][j] : 0;
		}

	dir == -1 ? max_s /= 255 : max_s = 1.0f;
	max_s = 1.0f / max_s;
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			OutputImage[i][j] *= max_s;

	if (dir == -1)
	{
		for (int i = M - 1; i >= M / 2; i--)
		{
			for (int j = N - 1; j >= 0; j--)
			{
				double temp = OutputImage[i][j];
				OutputImage[i][j] = OutputImage[M - i][j];
				OutputImage[M - i][j] = temp;
			}
		}
		for (int i = M - 1; i >= 0; i--)
		{
			for (int j = N - 1; j >= N / 2; j--)
			{
				double temp = OutputImage[i][j];
				OutputImage[i][j] = OutputImage[i][N - j];
				OutputImage[i][N - j] = temp;
			}
		}

	}

	for (int i = 0; i < M; ++i)
		free(out[i]);
	free(out);
}


void FT::FFT(double * pFreqReal, double * pFreqImag, std::complex<double> * x, long int row, long int n, long int start)
{
	const double PI = 3.1415926;
	int N = row;
	for (int i = 1, j = 0; i < N; ++i)
	{
		for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
		if (i > j) swap(&x[i], &x[j]);
	}

	/* dynamic programming */
	for (int k = 2; k <= N; k <<= 1)
	{
		double omega = -2.0 * PI / (1.0*k);
		//利用Eular's equation計算傅立葉之實虛數部分
		//InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
		//InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);

		std::complex<double> dtheta(cos(omega), sin(omega));
		// 每k個做一次FFT
		for (int j = 0; j < N; j += k)
		{
			// 前k/2個與後k/2的三角函數值恰好對稱，
			// 因此兩兩對稱的一起做。
			std::complex<double> theta(1, 0);
			for (int i = j; i < j + k / 2; i++)
			{
				std::complex<double> a = x[i];
				std::complex<double> b = x[i + k / 2] * theta;
				x[i] = a + b;
				x[i + k / 2] = a - b;
				theta *= dtheta;
			}
		}
	}
	for (int i = 0; i < N; ++i)
	{
		if (start == 1)
		{
			//x[i] = std::conj(x[i]);
			x[i] /= N;
		}
	}
}


void FT::LowpassFilter(double** Real, double** Img, double** filter)
{
}

void FT::HighpassFilter(double** Real, double** Img, double** filter)
{
}
