#pragma once

#define USE_THETA
#define SHOW_LIGHT_CAMERA


#include <iostream>

#include <iomanip>
#include <complex>
#include <cstdio>
#include <stdlib.h>
#include <fstream>
#include <time.h>
// Asynchronous
#include <thread>
#include <future>

#include "openCVlibs.h"

#ifdef USE_THETA
#include "ThetaWifiStream.h"
#else
// Input video from OpenCV
#endif // USE_THETA



#define M_PI       3.14159265358979323846   // pi from math.h

using namespace std;
using namespace cv;

uchar* DATA;
const int ch = 3; //CHANNELS
int const NThreads = std::thread::hardware_concurrency();
const int maxl = 2;
const int maxtheta = 1000;
const int maxphi = 1000;
const int simpsonflag = 0;
float theta[maxtheta];
float costheta[maxtheta];
float sintheta[maxtheta];
float sintheta2[maxtheta];
float sinthetafull[maxtheta];
complex<float> expphi[maxphi];
complex<float> expphistar[maxphi];
complex<float> expphim[maxphi][2 * maxl + 1];
float cosphi[maxphi];
float sinphi[maxphi];
float phi[maxphi];
float plm[maxl + 1][maxl + 1][maxtheta];

const float rt2 = 1.4142136;

const int _size = 1000;
int _size_x = _size; // number of rows
int _size_y = _size; // number of columns
const int maxsize = maxtheta;

//int mymaxl = maxl;
float floatfile[ch][maxsize][maxsize];
float lightcoeffs[ch][maxl + 1][2 * maxl + 1];
float img[ch][maxsize][maxsize];

float attenuation[maxl + 1];


void setuplegendre(void);
int Spherical_Harmonics();


bool reload_image = true; // A flag that signals reloading input image
string input_image_filename = "";

bool save_frame = false;
string save_frame_path = "";

float plmvalreal(int l, int m, float i)
{
	//assert(l >= 0 && l <= maxl && m >= -l && m <= l);
	int mindex = m; if (m < 0) mindex *= -1;
	int iint = (int)i;
	float fraci = i - iint;
	float t;
	//assert(i >= 0 && i <= maxtheta - 1);
	if (fabs(fraci)<1.0e-4)		t = plm[l][mindex][iint];
	else						t = fraci*plm[l][mindex][iint + 1] + (1 - fraci)*plm[l][mindex][iint];
	if (m != 0) t *= rt2;
	return t;
}

float combfloat(unsigned long n, unsigned long m)
{
	//assert(n >= 0 && m >= 0 && n >= m);
	if (m > n / 2) m = n - m;
	unsigned long num = 1, den = 1;
	int i;
	float retval = num / den;
	for (i = 1; i <= m; i++) {
		num = n - i + 1;
		den = i;
		retval *= (float)num / (float)den;
	}
	return retval;
}

float fact(float m)
{
	// Computes sqrt[2m!/(2^m m!)^2] 
	float retval = sqrt((2 * m + 1) / (4 * M_PI));
	for (int i = 1; i <= 2 * m - 1; i += 2)
		retval *= -sqrt((float)i / (float)(i + 1));
	return retval;
}

float expreal(int m, float j)
{
	int iint = (int)j;
	float fraci = j - iint;
	float t;
	//  assert(j >= 0 && j <= maxphi - 1) ;
	//assert(j >= 0 && j <= maxphi);
	if (fabs(fraci)<1.0e-4)
	{
		if (m >= 0) t = expphim[iint][m + maxl].real();
		else		t = expphim[iint][-m + maxl].imag();
	}
	else {
		if (m >= 0) t = (1 - fraci)*expphim[iint][m + maxl].real() + fraci *expphim[(iint + 1) % maxphi][m + maxl].real();
		else		t = (1 - fraci)*expphim[iint][-m + maxl].imag() + fraci *expphim[(iint + 1) % maxphi][-m + maxl].imag();
	}
	return t;
}

void setuplegendre(void)
{
	//  simpsonflag = 1 ;
	int i, j, p, q, k, l, m, n, s;
	for (i = 0; i < maxtheta; i++) {
		theta[i] = M_PI*i / maxtheta;
		costheta[i] = cos(theta[i]);
		sintheta[i] = sin(theta[i]);
		sinthetafull[i] = sin(theta[i] / 2.0);
		sintheta2[i] = sin(theta[i] / 4.0);
	}
	for (i = 0; i < maxphi; i++)
	{
		phi[i] = 2.0*M_PI*i / (maxphi);
		cosphi[i] = cos(phi[i]);
		sinphi[i] = sin(phi[i]);
		complex<float> z(cosphi[i], sinphi[i]);
		complex<float> zstar(cosphi[i], -sinphi[i]);
		expphi[i] = z;
		expphistar[i] = zstar;
		for (l = -maxl; l <= maxl; l++)
			expphim[i][l + maxl] = complex<float>(cos(l*phi[i]), sin(l*phi[i]));
	}
	for (m = 0; m <= maxl; m++) 
	{
		float fac = fact(m);
		for (i = 0; i < maxtheta; i++) 
		{
			plm[m][m][i] = fac*pow(sintheta[i], m);
			if (m < maxl) plm[m + 1][m][i] = sqrt(2 * m + 3) / (2 * m + 1)*costheta[i] * (2 * m + 1)*plm[m][m][i];
		}
		for (l = m + 2; l <= maxl; l++) 
		{
			for (i = 0; i < maxtheta; i++) 
			{
				float v1 = sqrt((2 * l + 1.0)*(2 * l - 1.0)*(l - m) / ((float)(l + m))) * costheta[i];
				float v2 = sqrt((2 * l + 1.0) / (2 * l - 3.0)*(l - m) / ((float)(l + m))*(l - m - 1.0)*(l + m - 1.0));
				plm[l][m][i] = (v1*plm[l - 1][m][i] - v2*plm[l - 2][m][i]) / (l - m);
			}
		}
	}
}

float readfloat(FILE * fp)
{
	float val;
	unsigned char *c = (unsigned char *)(&val);
	fscanf(fp, "%c%c%c%c", &(c[3]), &(c[2]), &(c[1]), &(c[0]));
	return val;
}

/*
void writeimage(char *filename, float scalefac)
{
	FILE *fp;
	int i, j, k;
	//assert(fp = fopen(filename, "wb"));
	fprintf(fp, "P6\n%d %d\n%d\n", _size_y, _size_x, 255);
	for (i = 0; i < _size_x; i++)
		for (j = 0; j < _size_y; j++)
			for (k = 0; k < ch; k++)
				fprintf(fp, "%c", DATA[ch*(_size_y*i + j) + 2 - k]);
	fclose(fp);
}

void writeimage(char *filename, float scalefac, float image[ch][maxsize][maxsize])
{
	FILE *fp;
	int i, j, k;
	//assert(fp = fopen(filename, "wb"));
	fprintf(fp, "P6\n%d %d\n%d\n", _size_y, _size_x, 255);
	for (i = 0; i < _size_x; i++)
		for (j = 0; j < _size_y; j++)
			for (k = 0; k < 3; k++)
			{
				float v = scalefac*image[k][i][j];
				unsigned char c;
				if (v < 0) c = (unsigned char)0;
				else if (v > 255) c = (unsigned char)255;
				else c = (unsigned char)v;
				fprintf(fp, "%c", c);
			}
	fclose(fp);
}

*/

float arrayval_phi(int channel, int ii, float index)
{
	int a, b;
	a = (int)index;
	b = a + 1;
	b = b % _size_y;

	//assert(a >= 0 && a < _size_y && b >= 0 && b < _size_y);
	float frac = index - a;
	return frac * (float)DATA[ch*(_size_y*ii + b) + 2 - channel] + (1 - frac)* (float)DATA[ch*(_size_y*ii + a) + 2 - channel];
}

float arrayval_theta(float fn[maxsize], float index)
{
	int a, b;
	a = (int)index; b = a + 1;
	if (b == _size_x) b = _size_x - 1;
	//assert(a >= 0 && a < _size_x && b >= 0 && b < _size_x);
	float frac = index - a;
	return frac * fn[b] + (1 - frac)*fn[a];
}

float integratephi(int channel, int m, int ii)
{

	float retval = 0;
	float mulfac = 2.0*M_PI / maxphi;

	for (int i = 0; i < maxphi; i++)
	{
		float iposn = (float)i * (float)_size_y / (float)maxphi;
		float simpsons = 1.0;
		if (simpsonflag)
		{
			if (i == 0) simpsons = 1.0 / 3.0;
			else if (i == maxphi - 1) simpsons = 1.0;
			else if (i == maxphi - 2) simpsons = 4.0 / 3.0;
			else if (i % 2 == 0) simpsons = 2.0 / 3.0;
			else simpsons = 4.0 / 3.0;
		}
		retval += expreal(m, i)*mulfac*arrayval_phi(channel, ii, iposn)*simpsons;
	}

	return retval;
}

float integratetheta(int l, int m, float fn[maxsize])
{
	float retval = 0;
	float mulfac = M_PI / maxtheta;
	float iposn, simpsons;
	for (int i = 0; i < maxtheta; i++)
	{
		iposn = (float)i * (float)_size_x / (float)maxtheta;
		simpsons = 1.0;
		if (simpsonflag)
		{
			if (i == 0) simpsons = 1.0 / 3.0;
			else if (i == maxtheta - 1) simpsons = 1.0;
			else if (i == maxtheta - 2) simpsons = 4.0 / 3.0;
			else if (i % 2 == 0) simpsons = 2.0 / 3.0;
			else simpsons = 4.0 / 3.0;
		}
		retval += arrayval_theta(fn, iposn)*sintheta[i] * mulfac*plmvalreal(l, m, i)*simpsons;
	}
	return retval;
}

void findcoeffslm(int channel)
{
	float fnt[2 * maxl + 1][maxtheta];
	int l, m, i;
	for (i = 0; i < _size_x; i++)
		for (m = -maxl; m <= maxl; m++)
			fnt[m + maxl][i] = integratephi(channel, m, i);
	for (l = 0; l <= maxl; l++)
		for (m = -l; m <= l; m++)
			lightcoeffs[channel][l][m + maxl] = integratetheta(l, m, fnt[m + maxl]);
}

void findcoeffslm()
{
	float fnt[ch][2 * maxl + 1][maxtheta];
	int l, m, i;

	for (i = 0; i < _size_x; i++)
		for (m = -maxl; m <= maxl; m++)
			for (int c = 0; c < 1; c++)
				fnt[c][m + maxl][i] = integratephi(c, m, i);
	for (int c = 0; c < ch; c++)
		for (l = 0; l <= maxl; l++)
			for (m = -l; m <= l; m++)
				lightcoeffs[c][l][m + maxl] = integratetheta(l, m, fnt[c][m + maxl]);
}

void findcoeffslm_async()
{
	float SCALE = 0.03;
	float fnt[ch][2 * maxl + 1][maxtheta];
	int l;
	int PerThread = _size_x / NThreads;
	vector<future<void>> futures;
	for (unsigned int iTh = 0; iTh < NThreads; ++iTh)
	{
		future<void> result(async([iTh, &PerThread, &fnt]()
		{
			int end = (iTh + 1) * PerThread;
			if (iTh == NThreads - 1) end = _size_x;
			for (int i = iTh * PerThread; i < end; i++)
				for (int channel = 0; channel < ch; channel++)
					for (int m = -maxl; m <= maxl; m++)
						fnt[channel][m + maxl][i] = integratephi(channel, m, i);
		}
		));

		futures.push_back(move(result));
	}
	for (vector<future<void>>::iterator it = futures.begin(); it < futures.end(); it++)
		it->get();
	for (int channel = 0; channel < ch; channel++)
		for (l = 0; l <= maxl; l++)
			for (int m = -l; m <= l; m++)
				lightcoeffs[channel][l][m + maxl] = SCALE * integratetheta(l, m, fnt[channel][m + maxl]);
}

void unfourier(float coeffs[maxl + 1][2 * maxl + 1], float vals[maxsize][maxsize])
{
	int i, j, l, m;

	float fnt[2 * maxl + 1][maxsize];

	for (i = 0; i <= 2 * maxl; i++)
		for (j = 0; j < _size_x; j++)
			fnt[i][j] = 0;

	for (m = -maxl; m <= maxl; m++)
	{
		int mmod = (m > 0) ? m : -m;
		for (l = mmod; l <= maxl; l++)
		{
			for (i = 0; i < _size_x; i++)
			{
				float iposn = (float)i* (float)maxtheta / (float)_size_x;
				if (iposn > maxtheta - 1) iposn = maxtheta - 1.0000001;
				fnt[m + maxl][i] += coeffs[l][m + maxl] * plmvalreal(l, m, iposn);
			}
		}
	}

	for (i = 0; i < _size_x; i++)
		for (j = 0; j < _size_y; j++)
		{
			vals[i][j] = 0;
			for (m = -maxl; m <= maxl; m++)
			{
				float jposn = (float)j* (float)maxphi / (float)_size_y;
				vals[i][j] += fnt[m + maxl][i] * expreal(m, jposn);
			}
		}

}

int Spherical_Harmonics()
{
	if (!DATA) return 0;

	int channel, l, m;

	findcoeffslm_async();
	return 0;
}

int webcam_spherical_harmonics()
{
#ifdef USE_THETA
	string thetaUrl = "http://192.168.1.1:80";
	RicohTheta theta(thetaUrl);
	//string SessionID = theta._startSession();
	//theta._getLivePreview(SessionID);
	theta.startSession();
	theta.startLivePreview();

#else
	VideoCapture cap(0); // open the default camera

	if (!cap.isOpened())  // check if we succeeded
		return -1;
#endif

	setuplegendre();

	Mat frame;
	int fPS, max_FPS, min_FPS;
	max_FPS = 0;
	min_FPS = 100000;
#ifdef SHOW_LIGHT_CAMERA
	namedWindow("Light Cam");
#endif
	for (;;)
	{
#ifdef USE_THETA
		theta.updateLiveFrame();
		frame = theta.getLiveFrame();
#else
		int start_s = clock();
		cap >> frame; // get a newframe from camera
#endif
		if (save_frame)
		{
			imwrite("web/"+save_frame_path, frame); // A JPG FILE IS BEING SAVED
			save_frame = false;
		}
#ifdef SHOW_LIGHT_CAMERA
		if (frame.data)
		{
			imshow("Light Cam", frame);
			waitKey(10);
		}
		else
		{
			printf("ERROR: No Frame Data To Show.\r");
		}

#endif

		resize(frame, frame, cv::Size::Size_(64,
#ifdef USE_THETA 
			32
#else 
			48
#endif
		), 0, 0, INTER_CUBIC);

		DATA = frame.data;
		_size_x = frame.rows;
		_size_y = frame.cols;
		Spherical_Harmonics();
		/*
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L00 ", lightcoeffs[0][0][2], lightcoeffs[1][0][2], lightcoeffs[2][0][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L1m1", lightcoeffs[0][1][1], lightcoeffs[1][1][1], lightcoeffs[2][1][1]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L10 ", lightcoeffs[0][1][2], lightcoeffs[1][1][2], lightcoeffs[2][1][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L11 ", lightcoeffs[0][1][3], lightcoeffs[1][1][3], lightcoeffs[2][1][3]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m2", lightcoeffs[0][2][0], lightcoeffs[1][2][0], lightcoeffs[2][2][0]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m1", lightcoeffs[0][2][1], lightcoeffs[1][2][1], lightcoeffs[2][2][1]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L20 ", lightcoeffs[0][2][2], lightcoeffs[1][2][2], lightcoeffs[2][2][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L21 ", lightcoeffs[0][2][3], lightcoeffs[1][2][3], lightcoeffs[2][2][3]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L22 ", lightcoeffs[0][2][4], lightcoeffs[1][2][4], lightcoeffs[2][2][4]);
		*/
	}

	// the camera will be deinitialized automatically in VideoCapture destructor
	return 0;
}


int file_spherical_harmonics()
{
	/*
	FILE * file = fopen(filename.c_str(), "r+");
	if (file == NULL) return -1;
	printf("File found.\n");
	fseek(file, 0, SEEK_END);
	long int size = THETA_FRAME_SIZE;// ftell(file);
	printf("File size:%i\n",size);
	fclose(file);
	// Reading data to array of unsigned chars
	file = fopen(filename.c_str(), "r+");
	uint8_t * Image = (uint8_t *)malloc(size);
	int bytes_read = fread(Image, sizeof(unsigned char), size, file);
	fclose(file);
	
	Mat rawData = Mat(1, THETA_FRAME_SIZE, CV_8UC1, Image);
	
	Mat img = imdecode(rawData, IMREAD_ANYCOLOR);
	*/
	Mat img = imread(input_image_filename, 1);
	setuplegendre();

	Mat frame;
	int fPS, max_FPS, min_FPS;
	max_FPS = 0;
	min_FPS = 100000;
	frame = img;
	int jpg_start_indx = 0;
	reload_image = false;
	while (!reload_image)
	{
		resize(frame, frame, cv::Size::Size_(64,
#ifdef USE_THETA 
			32
#else 
			48
#endif
		), 0, 0, INTER_CUBIC);

		DATA = frame.data;
		_size_x = frame.rows;
		_size_y = frame.cols;
		Spherical_Harmonics();
		/*
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L00 ", lightcoeffs[0][0][2], lightcoeffs[1][0][2], lightcoeffs[2][0][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L1m1", lightcoeffs[0][1][1], lightcoeffs[1][1][1], lightcoeffs[2][1][1]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L10 ", lightcoeffs[0][1][2], lightcoeffs[1][1][2], lightcoeffs[2][1][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L11 ", lightcoeffs[0][1][3], lightcoeffs[1][1][3], lightcoeffs[2][1][3]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m2", lightcoeffs[0][2][0], lightcoeffs[1][2][0], lightcoeffs[2][2][0]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L2m1", lightcoeffs[0][2][1], lightcoeffs[1][2][1], lightcoeffs[2][2][1]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L20 ", lightcoeffs[0][2][2], lightcoeffs[1][2][2], lightcoeffs[2][2][2]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L21 ", lightcoeffs[0][2][3], lightcoeffs[1][2][3], lightcoeffs[2][2][3]);
		printf("%s\t%10.5f\t%10.5f\t%10.5f\n", "L22 ", lightcoeffs[0][2][4], lightcoeffs[1][2][4], lightcoeffs[2][2][4]);
		*/
	}
	file_spherical_harmonics();
	// the camera will be deinitialized automatically in VideoCapture destructor
	return 0;
}
