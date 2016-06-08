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
	return fabs(P_f[s_idx] - P_f[t_idx]);
}

// For the 4-stencil discretization suggested in Bai and Sapiro 2008, and using
// the approximate gradient
// W_xy = |P_F(x) - P_F(y)|
// as the weight of edge (x, y), runs Dijkstra's algorithm to find the single-
// -source shortest paths from the source s.
double Dijkstra(double * P_f, int R, int C, Point * s, Point * t)
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
	/*for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			Point p;
			p.x = i;
			p.y = j;
			Q.push(make_pair(p, dist[getIndex(R, C, &p)]));
		}
	}*/
	Q.push(make_pair(*s, 0.0));

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
		Point * u = new Point;
		u->x = temp.x;
		u->y = temp.y;
		// If we are at the target, we can stop early.
		if ((u->x == t->x) && (u->y == t->y))
		{
			delete u;
			return dist[getIndex(R, C, t)];
		} 
		// EXTRACT-MIN frrom Q.
		Q.pop();
		//cerr << "This point has x coordinate " << u.x << " and y coordinate " << u.y << endl;
		// Find the adjacent Points to this Point.
		const int dx[4] = {-1, 0, 1, 0};
		const int dy[4] = {0, -1, 0, 1};
		for (int i = 0; i < 4; i++)
		{
			Point * v = new Point;
			v->x = u->x + dx[i];
			v->y = u->y + dy[i];
			//cerr << "This point has x coordinate " << v->x << " and y coordinate " << v->y << endl;
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
					//cerr << "About to Push!" << endl;
					Q.push(make_pair(*v, new_dist));
					//cerr << "Pushed!" << endl;
				}
			}
			delete v;
		}
		delete u;
	}

	//cerr << "Dijkstra's finished running!" << endl;

	return dist[getIndex(R, C, t)];
}

// For the 4-stencil discretization suggested in Bai and Sapiro 2008, and using
// the approximate gradient
// W_xy = |P_F(x) - P_F(y)|
// as the weight of edge (x, y), runs Dijkstra's algorithm to find the single-
// -source shortest paths from the source s.
double DijkstraR(double * P_f, int R, int C, Point * s, Point * t)
{
	//cerr << "Started running Dijkstra with souce (" << s->x << " , " << s->y << ")" << endl;
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
	// Set the distance of the source to 0.
	dist[s->x * C + s->y] = 0.0;

	// Now populate our priority queue.
	/*for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			Point p;
			p.x = i;
			p.y = j;
			Q.push(make_pair(p, dist[getIndex(R, C, &p)]));
		}
	}*/
	Q.push(make_pair(*s, 0.0));

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
		Point * u = new Point;
		u->x = temp.x;
		u->y = temp.y;
		// If we are at the target, we can stop early.
		if ((u->x == t->x) && (u->y == t->y))
		{
			delete u;
			return dist[t->x * C + t->y];
		} 
		// EXTRACT-MIN frrom Q.
		Q.pop();
		//cerr << "This point has x coordinate " << u.x << " and y coordinate " << u.y << endl;
		// Find the adjacent Points to this Point.
		const int dx[4] = {-1, 0, 1, 0};
		const int dy[4] = {0, -1, 0, 1};
		for (int i = 0; i < 4; i++)
		{
			Point * v = new Point;
			v->x = u->x + dx[i];
			v->y = u->y + dy[i];
			//cerr << "This point has x coordinate " << v->x << " and y coordinate " << v->y << endl;
			// Check to see if each neighbor is valid.
			if ((0 <= v->x && v->x < R) && (0 <= v->y && v->y < C))
			{
				//cerr << "Relax the edge." << endl;
				// Relax with respect to the edge (u, v).
				double w = fabs(P_f[u->x * C + u->y] - P_f[v->x * C + u->y]);
				double new_dist = dist[u->x * C + u->y] + w;
				if (new_dist < dist[v->x * C + v->y])
				{
					dist[v->x * C + v->y] = new_dist;
					// The <queue> library doesn't have a DECREASE-KEY
					// function, so we just add another QueueElem with the same
					// value but decreased key.
					//cerr << "About to Push!" << endl;
					Q.push(make_pair(*v, new_dist));
					//cerr << "Pushed!" << endl;
				}
			}
			delete v;
		}
		delete u;
	}

	//cerr << "Dijkstra's finished running!" << endl;

	return dist[t->x * C + t->y];
}

// Actually calculates the geodesic distance from a scribble to a given Point
// x. This is just the minimum over all points in the scribble 
double getDistance(double * P_f, int R, int C, vector<Point> scribble, Point x)
{
	// Run Dijkstra's algorithm with Point x as the source. For each Point in
	// the scribble, get the distance.
	//vector<double> dists = Dijkstra(P_f, R, C, &x);
	vector<double> sdists;
	for (int i = 0; i < scribble.size(); i++)
	{
		Point * temp = new Point;
		temp->x = scribble[i].x;
		temp->y = scribble[i].y;
		double dist = Dijkstra(P_f, R, C, temp, &x);
		sdists.push_back(dist);
		delete temp;
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

double getDistanceR(double * P_f, int R, int C, vector<Point> scribble, Point x)
{
	// Run Dijkstra's algorithm with Point x as the source. For each Point in
	// the scribble, get the distance.
	//vector<double> dists = Dijkstra(P_f, R, C, &x);
	vector<double> sdists;
	for (int i = 0; i < scribble.size(); i++)
	{
		Point * temp = new Point;
		temp->x = scribble[i].x;
		temp->y = scribble[i].y;
		double dist = DijkstraR(P_f, R, C, temp, &x);
		sdists.push_back(dist);
		delete temp;
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

	//cerr << "There are " << Q.size() << " elements in the queue" << endl;

	while(!Q.empty())
	{
		//err << "Starting main loop!" << endl;
		// Get the point at the top of the queue.
		Point u = Q.top().first;
		// Now pop the QueueElem.
		Q.pop();
		// Update all of this point's neighbors.
		const int dx[4] = {-1, 0, 1, 0};
		const int dy[4] = {0, -1, 0, 1};
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

// julienr's code FOR TEST PURPOSES
void GeodesicDistanceMap(const std::vector<Point>& sources,
double* height,
int W,
int H,
double* dists) {
// The algorithm is actually equivalent to running Dijkstra once for each
// source and then keeping the minimum distance.
// This is similar to "SHORTEST-PATH FOREST WITH TOPOLOGICAL ORDERING"
//
// But we do it all at once so it should be faster. It works as follow :
// 1. assign a distance of 0 to all source nodes, infinity to others
// 2. create a list of unvisited nodes consisting of the source nodes
// 3. visit neighbors of the current node
// - if the current node allows a shortest path to the neighbor
// - update the neighbor dist
// - add the neighbor to the unvisited nodes
// 4. remove current node from visited
// 5. pick the node with the smallest distance from the unvisited node as
// the new current
const int N = W*H;
// priority queue
typedef pair<int, double> PriorityEntry;
auto comp = [](const PriorityEntry& e1, const PriorityEntry& e2) {
return e1.second > e2.second;
};
priority_queue<PriorityEntry, vector<PriorityEntry>, decltype(comp)> Q(comp);
for (int i = 0; i < N; ++i) {
dists[i] = numeric_limits<double>::max();
}
for (const Point& p : sources) {
const int i = W*p.y + p.x;
dists[i] = 0;
Q.push(make_pair(i, 0));
}
//dx dy pairs for neighborhood exploration
const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = { 0, 1, 0, -1};
// main loop
while(!Q.empty()) {
int u = Q.top().first;
const int ux = u % W;
const int uy = u / W;
Q.pop();
// explore neighbors
for (int i = 0; i < 4; ++i) {
const int vx = ux + dx[i];
const int vy = uy + dy[i];
if ((vx < 0 || vx >= W) || (vy < 0 || vy >= H)) {
continue;
}
const int v = vy*W + vx;
const double w = fabs(height[v] - height[u]);
if ((dists[u] + w) < dists[v]) { // we found a shortest path to v
dists[v] = dists[u] + w;
// TODO: should UPDATE existing v (instead of duplicating)
Q.push(make_pair(v, dists[v]));
}
}
}
}


double * test(double * P_f, int R, int C, vector<Point> scribble)
{
	double * dists = new double[R * C];
	GeodesicDistanceMap(scribble, P_f, C, R, dists);
	return dists;
}