#include<bits/stdc++.h>

#include <bits/stdc++.h>

constexpr double EPS = 1e-7;


struct Point{
	Point(double arg_x,double arg_y)
        :x(arg_x), y(arg_y) {}

	Point() {x = y = 0.0;}

    Point operator + (const Point &p){ return Point(x+p.x,y+p.y); }
	Point operator - (const Point &p){ return Point(x-p.x,y-p.y);}
	Point operator * (const double &a){ return Point(a*x,a*y); }
	Point operator / (const double &a){ return Point(x/a,y/a); }
    Point& operator+= (const Point &p)
    {
        this->x += p.x;
        this->y += p.y;
        return *this;
    }

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

bool doubleCompare(double a, double b) 
{
    return fabs(a - b) <= EPS;
}


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
    vec2 base(xp2 - xp1, yp2 - yp1);
    int q;
    std::cin >> q;
    while (q--) {
        double xq, yq;
        std::cin >> xq >> yq;
        vec2 v0(xq - xp1, yq - yp1);
        double p = cross(base, v0);
        double len = abs(base) * abs(v0);
        double dotVal = dot(base, v0);
        // std::cout << std::fixed << std::setprecision(15);
        // std::cout << p << " " << len << " " << dotVal << "\n";
        if (doubleCompare(p, 0)) {
            double baseN = norm(base);
            double v0N = norm(v0);
            if (doubleCompare(dotVal, -len) && !doubleCompare(len, 0.0))
                std::cout << "ONLINE_BACK\n";
            else{
                if (v0N > baseN) 
                    std::cout << "ONLINE_FRONT\n";
                else
                    std::cout << "ON_SEGMENT\n";
            }
        }  else if (p > 0 ) {
            std::cout << "COUNTER_CLOCKWISE\n";
        } else if (p < 0) {
            std::cout << "CLOCKWISE\n";
        }
    }
    return 0;
}