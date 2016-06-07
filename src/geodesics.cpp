// Calculates the geodesic distances between the foreground and background
// scribbles in color space via a discrtization of our graph and Dijkstra's
// algorithm. We follow the paper Bai and Sapiro 2008.

#include "geodesics.h"
#include <cstdio>
#include <algorithm>
#include <queue>
#include <limits>
#include <vector>
#include <unordered_map>
#include <utility>

// Given an R x C array and an entry index (x, y), returns the index if the
// array is representedas a 1D row-major array.
int getIndex(int R, int C, Point * p)
{
	return p->x * C + p->y;
}

// Given two adjacent Points s and t, calculates the edge weight
// W_st = |P_F(s) - P_F(t)|
// between them. 
double getWeight(double * P_f, int R, int C, Point * s, Point * t)
{
	int s_idx = getIndex(R, C, s);
	int t_idx = getIndex(R, C, t);
	return abs(P_f[s_idx] - P_f[t_idx]);
}

// For the 4-stencil discretization suggested in Bai and Sapiro 2008, and using
// the approximate gradient
// W_xy = |P_F(x) - P_F(y)|
// as the weight of edge (x, y), runs Dijkstra's algorithm to find the single-
// -source shortest paths from the source s.
vector<double> Dijkstra(double * P_f, int R, int C, Point * s)
{
	//cerr << "Started running Dijkstra with souce (" << s->x << " , " << s->y << ")" << endl;
	// First, create our priority queue for use in Dijkstra's algorithm.
	typedef pair<Point, double> QueueElem;
	// The third coordinate serves as our value; so for convenience we will
	// write a comparison operator between QueueElems.
	auto comp = [](const QueueElem& q1, const QueueElem& q2) {
		return q1.second > q2.second;
	};
	// Now create our priority queue Q as a vector of QueueElems with
	// our QueueElem key comparison operator.
	priority_queue<QueueElem, vector<QueueElem>, decltype(comp)> Q(comp);

	// Now we can start the implementation of Dijkstra's algorithm. First we
	// maintain an array of distance estimates for each Point. We initialize
	// each distance to infinity.
	vector<double> dist (R * C, numeric_limits<double>::max());
	// Set the distance of the source to 0.
	dist[getIndex(R, C, s)] = 0.0;

	// Now populate our priority queue.
	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			Point p;
			p.x = i;
			p.y = j;
			Q.push(make_pair(p, dist[getIndex(R, C, &p)]));
		}
	}
	// Now for the main loop of Dijkstra's algorithm. While the queue is not
	// empty, we EXTRACT-MIN from the [min-]priority queue and then explore
	// the [at most] four neighbors p' of the Point p. We then look at the
	// distance estimate
	// dist[p] + weight(p, p')
	// If this distance is smaller than the current distance to p', we update
	// the distance array. We then push the updated
	while (!Q.empty())
	{
		// Get the Point at the top of the queue.
		Point temp = Q.top().first;
		Point * u;
		u->x = temp.x;
		u->y = temp.y;
		// EXTRACT-MIN frrom Q.
		Q.pop();
		//cerr << "This point has x coordinate " << u.x << " and y coordinate " << u.y << endl;
		// Find the adjacent Points to this Point.
		const int dx[4] = {-1, 0, 1, 0};
		const int dy[4] = {0, -1, 0, 1};
		for (int i = 0; i < 4; i++)
		{
			Point * v;
			v->x = u->x + dx[i];
			v->y = u->y + dy[i];
			//cerr << "This point has x coordinate " << v.x << " and y coordinate " << v.y << endl;
			// Check to see if each neighbor is valid.
			if ((0 <= v->x && v->x < R) && (0 <= v->y && v->y < C))
			{
				//cerr << "Relax the edge." << endl;
				// Relax with respect to the edge (u, v).
				double new_dist = dist[getIndex(R, C, u)] + getWeight(P_f, R, C, u, v);
				if (new_dist < dist[getIndex(R, C, v)])
				{
					dist[getIndex(R, C, v)] = new_dist;
					// The <queue> library doesn't have a DECREASE-KEY
					// function, so we just add another QueueElem with the same
					// value but decreased key.
					Q.push(make_pair(*v, new_dist));
				}
			}
		}
	}

	return dist;
}

// Actually calculates the geodesic distance from a scribble to a given Point
// x. This is just the minimum over all points in the scribble 
double getDistance(double * P_f, int R, int C, vector<Point> scribble, Point x)
{
	// Run Dijkstra's algorithm with Point x as the source. For each Point in
	// the scribble, get the distance.
	vector<double> dists = Dijkstra(P_f, R, C, &x);
	vector<double> sdists;
	for (int i = 0; i < scribble.size(); i++)
	{
		sdists.push_back(dists[getIndex(R, C, &scribble[i])]);
	}
	double min_dist = numeric_limits<double>::max();
	for (int i = 0; i < sdists.size(); i++)
	{
		if (sdists[i] < min_dist)
		{
			min_dist = sdists[i];
		}
	}
	return min_dist;
}

