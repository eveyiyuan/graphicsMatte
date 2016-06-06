#include "kde.h"
#include <cstdio>
#include <string>
#include <cmath>
#include <algorithm>

void calcKDE(const vector<double>& xs, const vector<double>& ys, 
			 const vector<double>& weights, double eps,
			 vector<double>* probs) {

	const int W = 1;
	const int M = ys.size();

	vector<double> norm_x(xs.size(), 0);
	for(size_t i = 0; i < xs.size(); i++) {
		norm_x[i] = (xs[i] - 128)/ 128.0;
	}
	const int d = 1;
	const int N = xs.size();
	const double h = 0.1;

	figtree(d, N, M, W, const_cast<double*>(norm_x.data()), h,
		    const_cast<double*>(weights.data()), const_cast<double*>(ys.data()),
		    eps, probs->data());
}

void bgforeProb(const double* p_fore, const double* p_bg, int W, int H, double* prob, bool bg) {
	for(int i = 0; i < W*H; i++) {
		if(p_fore[i] == 0 && p_bg[i] == 0) {
			prob[i] = 1;
		}
		else {
			if (bg) {
				prob[i] = p_bg[i]/ (p_fore[i] + p_bg[i]);
			}
			else {
				prob[i] = p_fore[i]/ (p_fore[i] + p_bg[i]);
			}
		}
	}
}


void generateProbs(double* probs, Mat input, const vector< vector<double> > xs, double eps) {
	// Calculate the PDFs (one per channel)
	vector<double> ys;
	for(int i = 0; i < 255; i++) {
		ys.push_back((i - 128)/128.0);
	}


	vector< vector<double> > pdfs(3);
	vector<double> chan1(255,0);
	pdfs[0] = chan1;
	vector<double> chan2(255, 0);
	pdfs[1] = chan2;
	vector<double> chan3(255, 0);
	pdfs[2] = chan3;

	// Get kernel density estimation
	for(int i = 0; i < 3; i++) {
		vector<double> weights(xs[i].size(), 1/(double)xs[i].size());
		calcKDE(xs[i], ys, weights, eps, &pdfs[i]);
	}

	// Get the probability of being in the foreground for each pixel in the image
	
	for(int r = 0; r < input.rows; r++) {
		for (int c = 0; c < input.cols; c++) {
			probs[input.cols*r + c] = pdfs[0][input.at<Vec3b>(r, c)[0]]*
									  pdfs[1][input.at<Vec3b>(r, c)[1]]*
									  pdfs[2][input.at<Vec3b>(r, c)[2]];
		}
	}
}