#include <opencv2/opencv.hpp>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>

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
void CallBackFunc(int event, int x, int y, int flags, void* userdata) {
	Mat* source = (Mat*) userdata;

	if (event == EVENT_LBUTTONDOWN) {
		cout << "Left button " << x << "," << y << endl;
		drawBg = true;
		// (*source).at<uchar>(y, x) = 0;
		// changeSquare(*source, x, y, 10, 0);
	}
	else if (event == EVENT_RBUTTONDOWN) {
		cout << "Right button " << x << "," << y << endl;
		drawFore = true;
		// (*source).at<uchar>(y, x) = 255;
		// changeSquare(*source, x, y, 10, 255);
	}
	else if (event == EVENT_MOUSEMOVE) {
		
		if (drawBg == true) {
			cout << "Mouse move " << x << "," << y << endl;
			(*source).at<uchar>(y, x) = 0;
			changeSquare(*source, x, y, 10, 0);
		}
		else if (drawFore == true) {
			(*source).at<uchar>(y, x) = 255;
			changeSquare(*source, x, y, 10, 255);			
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
	float thresh = atof(argv[2]);
	Mat input = imread(argv[1], CV_LOAD_IMAGE_COLOR);
	Mat grey = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);


	
	namedWindow("Display", WINDOW_AUTOSIZE);
	setMouseCallback("Display", CallBackFunc, &grey);
	while (1) {
		imshow("Display", grey);
		if(waitKey(1) != -1) {
			break;
		}
	}
	// imshow("Display", input);
	// waitKey(0);
	simpleMatte(Vec3b(200,240,0), input, grey, thresh);
	imshow("Display", grey);
	waitKey(0);
	return 0;
}