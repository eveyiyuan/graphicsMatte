#ifndef _GEODESIC_H_
#define _GEODESIC_H_

#include "kde.h"
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;


vector<double> getDists(double * P_f, int R, int C, vector<Point> scribble);


#endif //_GEODESIC_H_