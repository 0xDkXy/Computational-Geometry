#include <bits/stdc++.h>

constexpr double EPS = 0.000000001;


struct Point{
	Point(double arg_x,double arg_y)
        :x(arg_x), y(arg_y) {}

	Point() {x = y = 0.0;}

    Point operator + (const Point &p){ return Point(x+p.x,y+p.y); }
	Point operator - (const Point &p){ return Point(x-p.x,y-p.y);}
	Point operator * (const double &a){ return Point(a*x,a*y); }
	Point operator / (const double &a){ return Point(x/a,y/a); }

    friend std::ostream& operator<<(std::ostream &out, const Point &a)
    {
        out << std::fixed << std::setprecision(10) << a.x << " " << a.y;
        return out;
    }

	double abs(){ return sqrt(norm()); }
	double norm(){ return x*x + y*y; }

	bool operator<(const Point &p) const{
		return x != p.x? x < p.x: y < p.y;
	}

	bool operator == (const Point &p) const{
		return fabs(x-p.x) < EPS && fabs(y-p.y) < EPS;
	}

   	/*
	void debug(){

		printf("(%.3lf,%.3lf)\n",x,y);
	}*/

	double x,y;
};

struct Line{
	Line()=default;
	Line(Point a,Point b):p({a,b}){}
    Line(double x1, double y1, double x2, double y2): p({Point(x1,y1),Point(x2,y2)}) {}
	/*void outPut(){
		printf("(%.3lf,%.3lf)-(%.3lf,%.3lf)\n",p[0].x,p[0].y,p[1].x,p[1].y);
	}*/
	Point p[2];
};


using vec2 = Point;

double norm(vec2 &a) {return a.x*a.x + a.y*a.y;}

double abs(vec2 &a) {return sqrt(norm(a));}

double dot(vec2 &a, vec2 &b) {return a.x*b.x + a.y*b.y;}
double dot(vec2 &&a, vec2 &b) {return a.x*b.x + a.y*b.y;}

double cross(vec2 &a, vec2 &b) {return a.x*b.y - a.y*b.x;}

// 求投影坐标
Point projection(Line &l, Point &p)
{
    vec2 base = l.p[1] - l.p[0];
    double r = dot(p - l.p[0], base)/norm(base);
    return l.p[0]+base*r;
} 


signed main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    double xp1,yp1,xp2,yp2;
    std::cin >> xp1 >> yp1 >> xp2 >> yp2;
    Line base(xp1, yp1, xp2, yp2);
    int q;
    std::cin >> q;
    while (q--) {
        double xq, yq;
        std::cin >> xq >> yq;
        Point tmp(xq, yq);
        Point aim = projection(base, tmp);
        // std::cout << std::fixed << std::setprecision(10) << aim.x << " " << aim.y << "\n";
        std::cout << aim << "\n";
    }
    return 0;
}