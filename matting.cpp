#include <opencv2/opencv.hpp>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include "kde.h"

using namespace std;
using namespace cv;

// Calculates a numerical difference between 2 pixels
// represented by a vector of 3 color values
static float colorDiff(const Vec3b &pix1, const Vec3b &pix2) {
	float diff = 0;
	for(int i = 0; i < 3; i++) {
		diff += (pix1[i] - pix2[i]) ^ 2;
	}

	return sqrt(diff);
}
void changeSquare(Mat &input, int x, int y, int l, int val) {
	for(int r = x - l/2; r < x + l/2; r++) {
		for (int c = y - l/2; c < y + 1/2; c++) {
			input.at<uchar>(c, r) = val;
		}
	}

}
void simpleMatte(const Vec3b &bg, const Mat &input, Mat &output, float thresh) {
	int c = input.cols;
	int r = input.rows;

	for (int x = 0; x < c; x++) {
		for (int y = 0; y < r; y++) {
			const Vec3b color = input.at<Vec3b>(y,x);
			if(colorDiff(bg, color) < thresh) {
				output.at<uchar>(y,x) = 0;
			}
			else {
				output.at<uchar>(y,x) = 255;
			}
		}
	}
}

static bool drawBg = false;
static bool drawFore = false;
static vector<Point> bg;
static vector<Point> fore;
static vector< vector<double> > fore_data(3);
static vector< vector<double> > bg_data(3);
void CallBackFunc(int event, int x, int y, int flags, void* userdata) {
	Mat* source = (Mat*) userdata;

	if (event == EVENT_LBUTTONDOWN) {
		cout << "Left button marking bg " << x << "," << y << endl;
		bg_data[0].push_back((double)(*source).at<Vec3b>(y,x).val[0]);
		bg_data[1].push_back((double)(*source).at<Vec3b>(y,x).val[1]);
		bg_data[2].push_back((double)(*source).at<Vec3b>(y,x).val[2]);
		(*source).at<Vec3b>(y, x).val[0] = 0;
		(*source).at<Vec3b>(y, x).val[1] = 0;
		(*source).at<Vec3b>(y, x).val[2] = 0;
		bg.push_back(Point(x,y));
		drawBg = true;
	}
	else if (event == EVENT_RBUTTONDOWN) {
		cout << "Right button marking fore" << x << "," << y << endl;
		fore_data[0].push_back((double)(*source).at<Vec3b>(y,x).val[0]);
		fore_data[1].push_back((double)(*source).at<Vec3b>(y,x).val[1]);
		fore_data[2].push_back((double)(*source).at<Vec3b>(y,x).val[2]);
		(*source).at<Vec3b>(y, x).val[0] = 255;
		(*source).at<Vec3b>(y, x).val[1] = 255;
		(*source).at<Vec3b>(y, x).val[2] = 255;
		fore.push_back(Point(x,y));
		drawFore = true;
	}
	else if (event == EVENT_MOUSEMOVE) {
		
		if (drawBg == true) {
			bg_data[0].push_back((double)(*source).at<Vec3b>(y,x).val[0]);
			bg_data[1].push_back((double)(*source).at<Vec3b>(y,x).val[1]);
			bg_data[2].push_back((double)(*source).at<Vec3b>(y,x).val[2]);
			circle(*source, Point(x,y), 5, Scalar(0,0,0), -1);
			bg.push_back(Point(x,y));
		}
		else if (drawFore == true) {
			fore_data[0].push_back((double)(*source).at<Vec3b>(y,x).val[0]);
			fore_data[1].push_back((double)(*source).at<Vec3b>(y,x).val[1]);
			fore_data[2].push_back((double)(*source).at<Vec3b>(y,x).val[2]);
			circle(*source, Point(x,y), 5, Scalar(255,255,255), -1);
			fore.push_back(Point(x,y));			
		}
	}
	else if (event == EVENT_LBUTTONUP) {
		drawBg = false;
	}
	else if (event == EVENT_RBUTTONUP) {
		drawFore = false;
	}
}


int main(int argc, char**argv)
{
	Mat input = imread(argv[1], CV_LOAD_IMAGE_COLOR);
	Mat grey = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	Mat input2 = imread(argv[1], CV_LOAD_IMAGE_COLOR);


	
	namedWindow("Display", WINDOW_AUTOSIZE);
	setMouseCallback("Display", CallBackFunc, &input2);
	while (1) {
		imshow("Display", input2);
		if(waitKey(1) != -1) {
			break;
		}
	}

	double * fore_probs = new double[input.rows * input.cols];
	double * bg_probs = new double[input.rows * input.cols];

	generateProbs(fore_probs, input, fore_data, 0.001);
	generateProbs(bg_probs, input, bg_data, 0.001);

	double * P_Fx = new double[input.rows * input.cols];
	double * P_Bx = new double[input.rows * input.cols];

	bgforeProb(fore_probs, bg_probs, input.cols, input.rows, P_Fx, false);
	bgforeProb(fore_probs, bg_probs, input.cols, input.rows, P_Bx, true);

	for(unsigned int r = 0; r < input.rows; r++) {
		for(unsigned int c = 0; c < input.cols; c++) {
			if(P_Fx[input.cols*r + c] > P_Bx[input.cols*r + c]) {
				grey.at<uchar>(r, c) = 255;
			}
			else if (P_Fx[input.cols*r + c] < P_Bx[input.cols*r + c]){
				grey.at<uchar>(r, c) = 0;
			}
			else {
				grey.at<uchar>(r, c) = 255/2;
			}
		}
	}
	imshow("Display", grey);
	waitKey(0);
	
	return 0;
}