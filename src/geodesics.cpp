// Calculates the geodesic distances between the foreground and background
// scribbles in color space via a discrtization of our graph and Dijkstra's
// algorithm. We follow the paper Bai and Sapiro 2008.

#include "geodesics.h"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <queue>
#include <limits>
#include <vector>
#include <unordered_map>
#include <utility>



vector<double> getDists(double * P_f, int R, int C, vector<Point> scribble)
{
	// First, create our priority queue for use in Dijkstra's algorithm.
	typedef pair<Point, double> QueueElem;
	// The second coordinate serves as our value; so for convenience we will
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

	// Use each of our pixels in the scribble as sources. This allows us to
	// find the minimum number distance from ANY pixel to ANY source, allowing
	// us to get away with one run of Dijkstra's algorithm and get the geodesic
	// distance of any pixel.
	for (int i = 0; i < scribble.size(); i++)
	{
		// Set each source to be distance 0 - it is 0 distance away from itself, after all.
		//cerr << "Adding source!" << endl;
		dist[scribble[i].y * C + scribble[i].x] = 0.0;
		Q.push(make_pair(scribble[i], 0.0));
	}

	// The main loop.
	// While our priority queue is not empty, we EXTRACT-MIN from our priority
	// queue to get the element with the smallest distance. We then update the
	// distances of its neighbors by relaxing each edge and calling
	// DECREASE-KEY if we need to make an update.

	const int dx[4] = {-1, 0, 1, 0};
	const int dy[4] = {0, -1, 0, 1};

	while(!Q.empty())
	{
		//err << "Starting main loop!" << endl;
		// Get the point at the top of the queue.
		Point u = Q.top().first;
		// Now pop the QueueElem.
		Q.pop();
		// Update all of this point's neighbors.
		
		// Our graph is a 4-neighbor connected graph.
		for (int i = 0; i < 4; i++)
		{
			Point v;
			v.x = u.x + dx[i];
			v.y = u.y + dy[i];
			// Do bounds checking.
			if ((0 <= v.x && v.x < C) && (0 <= v.y && v.y < R))
			{
				//cerr << "Blah!" << endl;
				// Calculate our distance in our discretization.
				double w = fabs(P_f[u.y * C + u.x] - P_f[v.y * C + v.x]);
				//cerr << "The weight is " << w << endl;
				// Relax our edge.
				if (dist[u.y * C + u.x] + w < dist[v.y * C + v.x])
				{
					//cerr << "Updated!" << endl;
					dist[v.y * C + v.x] = dist[u.y * C + u.x] + w;
					// Since the <queue> library does not have a DECREASE-KEY
					// function, add a duplicate instead.
					Q.push(make_pair(v, dist[u.y * C + u.x] + w));
				}
			}
		}
	}
	// Return our vector of distances.
	return dist;
}