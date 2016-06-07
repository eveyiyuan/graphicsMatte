#ifndef _GEODESIC_H_
#define _GEODESIC_H_

#include "kde.h"
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// Helper function to transform to 1D indices.
int getIndex(int R, int C, Point * p);

// Gets the edge weight in our discretization.
double getWeight(double * P_f, int R, int C, Point * s, Point * t);

// Run's Dijkstra's single-source shortest path algorithm in our 4-stencil
// discretization of our image from source s.
double Dijkstra(double * P_f, int R, int C, Point * s, Point * t);

// Calculates the actual geodesic distance from a given pixel x to a scribble
// (either the foreground or the background scribble).
double getDistance(double * P_f, int R, int C, vector<Point> scribble, Point x);

#endif //_GEODESIC_H_