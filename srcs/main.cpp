#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
#include <ctime>

#include <cmath>
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "gtx/scalar_multiplication.hpp"

#define SDL_MAIN_HANDLED
#include "SDL2_gfxPrimitives.h"

#define SCR_SIZE 700
#define SCALE 1.0

#define POINT_RADIUS 2
#define LINE_THICKNESS 1
#define COLOR_BLACK 0xFF000000
#define COLOR_GRAY 0xFFD0D0D0
#define COLOR_WHITE 0xFFFFFFFF
#define COLOR_RED 0xFF0000FF

#define MIN_SIZE 2
#define GROWTH_COEF 5

inline float length (float x, float y)
{
	return sqrt(x*x + y*y);
}

inline int toScr (float coord)
{
	return (coord*SCALE+1.0) * SCR_SIZE/2;
}

class Point {
public:
	float x, y;

	void draw (SDL_Renderer *renderer, Uint32 color);

	Point operator= (Point &point);

	friend bool operator== (Point &a, Point &b);
	friend bool operator!= (Point &a, Point &b);
};

void Point::draw (SDL_Renderer *renderer, Uint32 color)
{
	filledCircleColor(renderer, 
					  toScr(x), 
					  toScr(-y), 
					  POINT_RADIUS, 
					  color);
	aacircleColor(renderer, 
				  toScr(x), 
				  toScr(-y), 
				  POINT_RADIUS, 
				  color);
}

Point Point::operator= (Point &point)
{
	x = point.x;
	y = point.y;

	return *this;
}

bool operator== (Point &a, Point &b)
{
	return (a.x == b.x && a.y == b.y);
}

bool operator!= (Point &a, Point &b)
{
	return (a.x != b.x || a.y != b.y);
}

class Triangle {
public:
	std::vector<Point> vertices;
	std::vector<Triangle *> neighbours;

	bool visited;

	Triangle (Point a, Point b, Point c);
	~Triangle ();

	Point getCenter ();
	Point getCircumcenter ();

	void addNeighbour (Triangle *t, Point rib1, Point rib2);
	void removeNeighbours ();

	void draw (SDL_Renderer *renderer, Uint32 color, bool debug);

	Triangle operator= (Triangle &t);
};

Triangle::Triangle (Point a, Point b, Point c) : visited(false)
{
	vertices.push_back(a);
	vertices.push_back(b);
	vertices.push_back(c);
	for (int i = 0; i < 3; i++)
		neighbours.push_back(nullptr);
}

Triangle::~Triangle () 
{
	removeNeighbours();
}

Point Triangle::getCenter ()
{
	Point center;

	center.x = 0.0;
	center.y = 0.0;
	for (int i = 0; i < 3; i++) {
		center.x += vertices[i].x;
		center.y += vertices[i].y;
	}
	center.x /= 3.0;
	center.y /= 3.0;

	return center;
}

Point Triangle::getCircumcenter ()
{
	Point center;

	Point A = vertices[0];
	Point B = vertices[1];
	Point C = vertices[2];

	float D = 2.0*(A.x*(B.y - C.y) + B.x*(C.y - A.y) + C.x*(A.y - B.y));
	center.x = ((A.x*A.x + A.y*A.y)*(B.y-C.y)
				+ (B.x*B.x + B.y*B.y)*(C.y-A.y)
				+ (C.x*C.x + C.y*C.y)*(A.y-B.y))/D;
	center.y = ((A.x*A.x + A.y*A.y)*(C.x-B.x)
				+ (B.x*B.x + B.y*B.y)*(A.x-C.x)
				+ (C.x*C.x + C.y*C.y)*(B.x-A.x))/D;

	return center;
}

void Triangle::addNeighbour (Triangle *t, Point rib1, Point rib2)
{
	int i;
	for (i = 0; i < 3; i++)
		if (vertices[i] != rib1 && vertices[i] != rib2)
			break;

	int j;
	for (j = 0; j < 3; j++)
		if (t->vertices[j] != rib1 && t->vertices[j] != rib2)
			break;

	neighbours[i] = t;
	t->neighbours[j] = this;
}

void Triangle::removeNeighbours ()
{
	for (int i = 0; i < 3; i++)
		if (neighbours[i] != nullptr) {
			for (int j = 0; j < 3; j++) 
				if (neighbours[i]->vertices[j] != vertices[(i+1)%3]
					&& neighbours[i]->vertices[j] != vertices[(i+2)%3]) {
					
					neighbours[i]->neighbours[j] = nullptr;
				}
		}
}


void Triangle::draw (SDL_Renderer *renderer, Uint32 color, bool debug)
{
	aalineColor(renderer, 
				toScr(vertices[0].x), 
				toScr(-vertices[0].y), 
				toScr(vertices[1].x), 
				toScr(-vertices[1].y),
				color);

	aalineColor(renderer, 
				toScr(vertices[1].x), 
				toScr(-vertices[1].y), 
				toScr(vertices[2].x), 
				toScr(-vertices[2].y),
				color);

	aalineColor(renderer, 
				toScr(vertices[2].x), 
				toScr(-vertices[2].y), 
				toScr(vertices[0].x), 
				toScr(-vertices[0].y),
				color);

	if (debug) {
		Point circumcenter = getCircumcenter();
		Point radius;
		radius.x = circumcenter.x - vertices[0].x;
		radius.y = circumcenter.y - vertices[0].y;
		aacircleColor(renderer, 
					  toScr(circumcenter.x), 
					  toScr(-circumcenter.y), 
					  length(radius.x, radius.y) * SCALE * SCR_SIZE/2, 
					  COLOR_GRAY);
	}
}

Triangle Triangle::operator= (Triangle &t)
{
	for (int i = 0; i < 3; i++) {
		vertices[i] = t.vertices[i];
		neighbours[i] = t.neighbours[i];
	}

	return *this;
}

class Triangulation {
private:
	std::vector<Point> superstruct;

	void clearVisits (Triangle *t);
	void deleteRec (Triangle *t);
	void resizeCache ();
	Triangle *find (Point p);
	void setCache(Triangle *t, Point p);
	Triangle *iterate (Triangle *t, Point p);
	bool test (Triangle *t1, Triangle *t2, Point rib1, Point rib2);
	void flip (Triangle *t1, Triangle *t2, Point rib1, Point rib2);
	void removeSuperstructureRec (Triangle *t);

	void drawRec (SDL_Renderer *renderer, Triangle *t, bool debug);

public:
	Triangle *start;

	Triangle **cache;
	int cacheSize;
	int nPoints;

	Triangulation ();
	~Triangulation ();

	void addPoint (Point p);
	void removeSuperstructure ();

	void draw (SDL_Renderer *renderer, bool debug);
};

void Triangulation::clearVisits (Triangle *t)
{
	if (t != nullptr && t->visited) {
		t->visited = false;
		for (int i = 0; i < 3; i++)
			clearVisits(t->neighbours[i]);
	}
}

void Triangulation::deleteRec (Triangle *t)
{
	if (t != nullptr && !t->visited) {
		t->visited = true;
		for (int i = 0; i < 3; i++)
			deleteRec(t->neighbours[i]);
		delete t;
	}
}

void Triangulation::resizeCache ()
{
	Triangle **cacheOld = cache;
	cache = new Triangle *[cacheSize * cacheSize * 4];
	for (int i = 0; i < cacheSize; i++)
		for (int j = 0; j < cacheSize; j++) {
			cache[i*2*cacheSize + j*2] = cacheOld[i * cacheSize + j];
			cache[(i*2+1)*cacheSize + j*2] = cacheOld[i * cacheSize + j];
			cache[i*2*cacheSize + j*2+1] = cacheOld[i * cacheSize + j];
			cache[(i*2+1)*cacheSize + j*2+1] = cacheOld[i * cacheSize + j];
		}
	cacheSize *= 2;
	delete[] cacheOld;
}

bool Triangulation::test (Triangle *t1, Triangle *t2, Point rib1, Point rib2)
{
	int i;
	for (i = 0; i < 3; i++)
		if (t1->vertices[i] != rib1 && t1->vertices[i] != rib2)
			break;
	Point B = t1->vertices[i];
	for (int k = 0; k < superstruct.size(); k++)
		if (B == superstruct[k])
			return true;

	int j;
	for (j = 0; j < 3; j++)
		if (t2->vertices[j] != rib1 && t2->vertices[j] != rib2)
			break;
	Point A = t2->vertices[j];
	for (int k = 0; k < superstruct.size(); k++)
		if (A == superstruct[k])
			return true;

	float cA = (A.x-rib1.x) * (A.x-rib2.x) + (A.y-rib1.y) * (A.y-rib2.y);
	float cB = (B.x-rib1.x) * (B.x-rib2.x) + (B.y-rib1.y) * (B.y-rib2.y);

	if (cA < 0 && cB < 0) {
		return false;
	} else if (cA >= 0 && cB >= 0) {
		return true;
	} else {

		float sA = abs((A.x-rib1.x)*(A.y-rib2.y) - (A.x-rib2.x)*(A.y-rib1.y));
		float sB = abs((B.x-rib2.x)*(B.y-rib1.y) - (B.x-rib1.x)*(B.y-rib2.y));

		if (sA * cB + cA * sB >= 0) {
			// make shure we have a convex triangulation
			Point S;
			bool isSuperstruct = false;
			for (int i = 0; i < superstruct.size(); i++) {
				if (rib1 == superstruct[i]) {
					S = rib1;
					isSuperstruct = true;
					break;
				} else if (rib2 == superstruct[i]) {
					S = rib2;
					isSuperstruct = true;
					break;
				}
			}

			if (isSuperstruct) {
				float lS = length(S.x-A.x, S.y-A.y)*length(S.x-B.x, S.y-B.y);
				float cS = ((S.x-A.x)*(S.x-B.x) + (S.y-A.y)*(S.y-B.y)) / lS;
				
				float lA = length(A.x-rib1.x, A.y-rib1.y)
						   * length(A.x-rib2.x, A.y-rib2.y);
				float lB = length(B.x-rib1.x, B.y-rib1.y)
						   * length(B.x-rib2.x, B.y-rib2.y);
				float sumAB = (cA * cB - sA * sB) / (lA * lB);
	
				if (sumAB < -cS)
					return false;
			}
			return true;
		} else {
			return false;
		}
	}
}

void Triangulation::flip (Triangle *t1, Triangle *t2, Point rib1, Point rib2)
{
	if (!test(t1, t2, rib1, rib2)) {
		int t1A, t1Rib1, t1Rib2;
		for (int i = 0; i < 3; i++) {
			if (t1->vertices[i] == rib1)
				t1Rib1 = i;
			else if (t1->vertices[i] == rib2)
				t1Rib2 = i;
			else
				t1A = i;
		}

		int t2B, t2Rib1, t2Rib2;
		for (int i = 0; i < 3; i++) {
			if (t2->vertices[i] == rib1)
				t2Rib1 = i;
			else if (t2->vertices[i] == rib2)
				t2Rib2 = i;
			else
				t2B = i;
		}

		// external neighbours
		Triangle *n1, *n2, *n3, *n4;
		n1 = t1->neighbours[t1Rib2];
		n2 = t2->neighbours[t2Rib2];
		n3 = t1->neighbours[t1Rib1];
		n4 = t2->neighbours[t2Rib1];

		// remove all neighbours
		t1->removeNeighbours();
		t2->removeNeighbours();

		// veritces of new diagonal
		Point a = t1->vertices[t1A];
		Point b = t2->vertices[t2B];

		// creane "new" triangles (flip diagonal)
		t1->vertices[0] = rib2;
		t1->vertices[1] = b;
		t1->vertices[2] = a;

		t2->vertices[0] = rib1;
		t2->vertices[1] = a;
		t2->vertices[2] = b;

		// add neighbours to "new" triangles
		t1->neighbours[0] = t2;
		t2->neighbours[0] = t1;

		if (n3 != nullptr)
			t1->addNeighbour(n3, rib2, a);
		if (n4 != nullptr)
			t1->addNeighbour(n4, b, rib2);
		if (n1 != nullptr)
			t2->addNeighbour(n1, a, rib1);
		if (n2 != nullptr)
			t2->addNeighbour(n2, rib1, b);

		// flip external neighbours if needed 
		// (dont flip neighbours of "old" t1 triangle 
		// because we assume they already tested 
		// or will be tested lately)
		if (n2 != nullptr)
			flip(t2, n2, rib1, b);
		if (n4 != nullptr)
			flip(t1, n4, b, rib2);
	} 
}

void Triangulation::removeSuperstructureRec (Triangle *t)
{
	if (t != nullptr && !t->visited) {
		t->visited = true;
		bool superstructPart = false;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < superstruct.size(); j++)
				if (t->vertices[i] == superstruct[j]) {
					superstructPart = true;
					break;
				}
			if (superstructPart)
				break;
		}
		// move start triangle if its a superstructure triangle
		if (superstructPart && t == start)
			for (int i = 0; i < 3; i++)
				if (t->neighbours[i] != nullptr) {
					start = t->neighbours[i];
					break;
				}

		for (int i = 0; i < 3; i++) {
			if (t->neighbours[i] != nullptr)
				removeSuperstructureRec(t->neighbours[i]);
		}

		if (superstructPart) 
			delete t;
	}
}

void Triangulation::drawRec (SDL_Renderer *renderer, Triangle *t, bool debug)
{
	if (t != nullptr && !t->visited) {
		t->visited = true;
		bool superstructPart = false;
		for (int i = 0; i < 3; i++) {
			if (t->neighbours[i] != nullptr)
				drawRec(renderer, t->neighbours[i], debug);
			if (debug)
				for (int j = 0; j < superstruct.size(); j++)
					if (t->vertices[i] == superstruct[j]) {
						superstructPart = true;
						break;
					}
		}
		if (!superstructPart)
			t->draw(renderer, COLOR_BLACK, debug);
	}
}

Triangle *Triangulation::find (Point p)
{
	int i = p.y / (float)cacheSize;
	int j = p.x / (float)cacheSize;

	return cache[i*cacheSize + j];
}

void Triangulation::setCache (Triangle *t, Point p)
{
	int i = p.y / (float)cacheSize;
	int j = p.x / (float)cacheSize;

	cache[i*cacheSize + j] = t;
}

Triangle *Triangulation::iterate (Triangle *t, Point p)
{
	for (int i = 0; i < 3; i++) {
		Point vertex = t->vertices[i];

		Point rib1 = t->vertices[(i+1)%3];
		Point rib2 = t->vertices[(i+2)%3];

		glm::vec3 s1 = {rib1.x-vertex.x, rib1.y-vertex.y, 0.0};
		glm::vec3 s2 = {rib1.x-p.x, rib1.y-p.y, 0.0};
		glm::vec3 line = {rib1.x-rib2.x, rib1.y-rib2.y, 0.0};
		glm::vec3 prod1 = glm::cross(line, s1);
		glm::vec3 prod2 = glm::cross(line, s2);
		if (prod1.z * prod2.z < 0)
			return iterate(t->neighbours[i], p);
	}

	return t;
}

Triangulation::Triangulation () : 
	start(nullptr), cacheSize(MIN_SIZE), nPoints(0)
{
	// triangle which contains 2 by 2 screen square 
	// in which all points are located
	superstruct.push_back({-2.15, -1.0});
	superstruct.push_back({2.15, -1.0});
	superstruct.push_back({0.0, 2.73});

	start = new Triangle(superstruct[0], superstruct[1], superstruct[2]);

	cache = new Triangle *[cacheSize * cacheSize];
	for (int i = 0; i < cacheSize; i++)
		for (int j = 0; j < cacheSize; j++)
			cache[i*cacheSize + j] = start;
}

Triangulation::~Triangulation ()
{
	delete[] cache;

	deleteRec(start);
}

void Triangulation::addPoint (Point p)
{
	Triangle *target = iterate(find(p), p);

	float l[3];
	for (int i = 0; i < 3; i++)
		l[i] = length(p.x-target->vertices[i].x, p.y-target->vertices[i].y);

	bool inside = true;
	int index;
	for (index = 0; index < 3; index++) {
		float sL = length(target->vertices[index].x
						  -target->vertices[(index+1)%3].x,
						  target->vertices[index].y
						  -target->vertices[(index+1)%3].y);
		if (sL >= l[index] + l[(index+1)%3]) {
			inside = false;
			break;
		}
	}

	if (inside) {
		// if p located inside triangle
		// points of old triangle
 		Point a = target->vertices[0];
		Point b = target->vertices[1];
		Point c = target->vertices[2];
	
		// external neighbours
		Triangle *n1, *n2, *n3;
		n1 = target->neighbours[0];
		n2 = target->neighbours[1];
		n3 = target->neighbours[2];
	
		//remove all neigbours
		target->removeNeighbours();
	
		// create new triangles
		target->vertices[2] = p;
		Triangle *t1 = new Triangle(b, c, p);
		Triangle *t2 = new Triangle(c, a, p);
	
		// add neighbours to new triangles
		target->neighbours[0] = t1;
		target->neighbours[1] = t2;
		t1->neighbours[0] = t2;
		t1->neighbours[1] = target;
		t2->neighbours[0] = target;
		t2->neighbours[1] = t1;
	
		if (n3 != nullptr)
			target->addNeighbour(n3, a, b);
		if (n1 != nullptr)
			t1->addNeighbour(n1, b, c);
		if (n2 != nullptr)
			t2->addNeighbour(n2, c, a);
	
		// flip external neigbours if needed
		if (n3 != nullptr)  
			flip(target, n3, a, b);
		if (n1 != nullptr)
			flip(t1, n1, b, c);
		if (n2 != nullptr)
			flip(t2, n2, c, a);
	} else {
		// if p located on one of the sides
		// points of old triangle
 		Point a = target->vertices[index];
		Point b = target->vertices[(index+1)%3];
		Point c = target->vertices[(index+2)%3];
		
		// external neighbours of old triangle
		Triangle *n1, *n2;
		n1 = target->neighbours[index];
		n2 = target->neighbours[(index+1)%3];

		// external triangle
		Triangle *ext = target->neighbours[(index+2)%3];

		// external point and neighbours of external triangle
		Point d;
		Triangle *n3, *n4;
		for (int i = 0; i < 3; i++)
			if (ext->vertices[i] == a)
				n3 = ext->neighbours[i];
			else if (ext->vertices[i] == b)
				n4 = ext->neighbours[i];
			else
				d = ext->vertices[i];
	
		//remove all neigbours
		target->removeNeighbours();
		ext->removeNeighbours();
	
		// create new triangles
		target->vertices[0] = a;
		target->vertices[1] = c;
		target->vertices[2] = p;
		ext->vertices[0] = a;
		ext->vertices[1] = d;
		ext->vertices[2] = p;
		Triangle *t1 = new Triangle(b, c, p);
		Triangle *t2 = new Triangle(b, d, p);
	
		// add neighbours to new triangles
		target->neighbours[0] = t1;
		target->neighbours[1] = ext;
		ext->neighbours[0] = t2;
		ext->neighbours[1] = target;
		t1->neighbours[0] = target;
		t1->neighbours[1] = t2;
		t2->neighbours[0] = ext;
		t2->neighbours[1] = t1;
	
		if (n1 != nullptr)
			t1->addNeighbour(n1, b, c);
		if (n2 != nullptr)
			target->addNeighbour(n2, c, a);
		if (n3 != nullptr)
			t2->addNeighbour(n3, d, b);
		if (n4 != nullptr)
			ext->addNeighbour(n4, a, d);
	
		// flip external neigbours if needed
		if (n1 != nullptr)
			flip(t1, n1, b, c);
		if (n2 != nullptr)
			flip(target, n2, c, a);
		if (n3 != nullptr)
			flip(t2, n3, d, b);
		if (n4 != nullptr)
			flip(ext, n4, a, d);
	}

	setCache(target, p);
	nPoints++;
	if (nPoints >= GROWTH_COEF * cacheSize * cacheSize)
		resizeCache();
}

void Triangulation::removeSuperstructure ()
{
	removeSuperstructureRec(start);
	clearVisits(start);
}

void Triangulation::draw (SDL_Renderer *renderer, bool debug)
{
	drawRec(renderer, start, debug);
	clearVisits(start);
}

bool readPoints (std::string path, std::vector<Point> &points)
{
	std::ifstream fin(path);
	
	if (fin.is_open()) {
		bool first = true;
		float max = 0.0;
		float min = 0.0;
		
		while (fin.good()) {
			float x, y;
			fin >> x >> y;
			if (!fin.good() && !fin.eof())
				return false;

			if (first) {
				first = false;
				max = std::max(x, y);
				min = std::min(x, y);
			} else {
				if (x > max)
					max = x;
				else if (x < min)
					min = x;
				if (y > max)
					max = y;
				else if (y < min)
					min = y;
			}

			points.push_back({x, y});
		}
		fin.close();

		float scale = abs(max) + abs(min);
		float center = 0.0 - ((max + min) / 2.0);
		for (int i = 0; i < points.size(); i++) {
			points[i].x = ((points[i].x + center) / scale) * 1.9;
			points[i].y = ((points[i].y + center) / scale) * 1.9;
		}

	} else {
		return false;
	}

	return true;
}

int main ()
{
	// initial settings
	std::cout << "enter a path to file containing points\n> ";
	std::vector<Point> points;
	std::string path;
	std::cin >> path;
	if (!readPoints(path, points)) {
		std::cout << "ERROR::unable to open file or data is incorrect\n";
		return 1;
	}
	if (points.size() < 3) {
		std::cout << "ERROR::not enough points to build triangulation\n";
		return 1;
	}

	bool debug = false;
	std::cout << "would you like to enable debug mode:\n"
			  << "(enables step by step triangulation building and "
			  << "outputting circumcircles of triangles) "
			  << "(y/n)?\n" << "> ";
	char action;
	std::cin >> action;
	if (std::tolower(action) == 'y')
		debug = true;

	// setup window
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0){
		std::cout << "ERROR::unable to initialise SDL: " 
				  << SDL_GetError() << std::endl;
		return 1;
	}

	SDL_Window *window = SDL_CreateWindow("visualization", 
										  0, 25, 
										  SCR_SIZE, SCR_SIZE, 
										  SDL_WINDOW_SHOWN);
	if (window == nullptr){
		std::cout << "ERROR::unable to create window: " 
				  << SDL_GetError() << '\n';
		return 1;
	}

	SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 
												SDL_RENDERER_ACCELERATED 
												| SDL_RENDERER_PRESENTVSYNC);
	if (renderer == nullptr){
		std::cout << "ERROR::unable to create renderer: " 
				  << SDL_GetError() << '\n';
		return 1;
	}
	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "1");

	// build triangulation
	Triangulation triangulation;
	
	bool done = false;
	
	if (debug) {
		bool next = false;
		int index = 0;
		while (!done && index < points.size()) {
			SDL_Event ev;
			while (SDL_PollEvent( &ev )) {
				if ((SDL_QUIT == ev.type) 
					|| (SDL_KEYDOWN == ev.type 
					&& SDL_SCANCODE_ESCAPE == ev.key.keysym.scancode)) {
	
					done = true;
					break;
				}
				if (SDL_KEYDOWN == ev.type 
					&& SDL_SCANCODE_RETURN == ev.key.keysym.scancode) {
						
					next = true;
					break;
				}
			}
	
			if (next) {
				next = false;
				Point p = points[index++];
				triangulation.addPoint(p);
	
				SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
				SDL_RenderClear(renderer);
				for (int i = 0; i < index; i++) {
					points[i].draw(renderer, COLOR_BLACK);
				}
				triangulation.draw(renderer, debug);
				SDL_RenderPresent(renderer);
			}
		}
		triangulation.removeSuperstructure();
	} else {
		unsigned int start = clock();
		for (int i = 0; i < points.size(); i++) {
			Point p = points[i];
			triangulation.addPoint(p);
		}
		triangulation.removeSuperstructure();
		unsigned int end = clock();
		std::cout << "triangulation built in "
				  << (float)(end - start) / (float)CLOCKS_PER_SEC
				  << " seconds\n";
	}

	if (!done) {
		std::cout << "triangulation built successfully\n";
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
		SDL_RenderClear(renderer);
		for (int i = 0; i < points.size(); i++) {
			points[i].draw(renderer, COLOR_BLACK);
		}
		triangulation.draw(renderer, debug);
		SDL_RenderPresent(renderer);
	
		while (!done) {
			SDL_Event ev;
			while (SDL_PollEvent( &ev )) 
				if ((SDL_QUIT == ev.type) 
					|| (SDL_KEYDOWN == ev.type 
					&& SDL_SCANCODE_ESCAPE == ev.key.keysym.scancode)) {
	
					done = true;
					break;
				}
		}
	}

	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}