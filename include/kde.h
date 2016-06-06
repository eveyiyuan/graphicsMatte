#ifndef _KDE_H_
#define _KDE_H_

#include "figtree.h"
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void calcKDE(const vector<double>& xs, const vector<double>& ys, 
			 const vector<double>& weights, double eps,
			 vector<double>* probs);

void bgforeProb(const double* p_fore, const double* p_bg, int W, int H, double* prob, bool bg);

void generateProbs(double* probs, Mat input, const vector< vector<double> > xs, double eps);

#endif //_KDE_H_