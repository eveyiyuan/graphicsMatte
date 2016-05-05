#include <opencv2/opencv.hpp>
#include <math.h>
#include <stdlib.h>

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

int main(int argc, char**argv)
{
	float thresh = atof(argv[2]);
	Mat input = imread(argv[1], CV_LOAD_IMAGE_COLOR);
	Mat output = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);

	simpleMatte(Vec3b(200,240,0), input, output, thresh);
	namedWindow("Display", WINDOW_AUTOSIZE);
	imshow("Display", output);
	waitKey(0);

	return 0;
}