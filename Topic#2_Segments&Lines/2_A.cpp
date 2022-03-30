#include<bits/stdc++.h>

namespace CompG{
        
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

    bool doubleEqual(double a, double b) 
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
};


int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    int q;
    std::cin >> q;
    while (q--) {
        double x1,x2,x3,x4,y1,y2,y3,y4;
        std::cin >> x1 >> y1;
        std::cin >> x2 >> y2;
        std::cin >> x3 >> y3;
        std::cin >> x4 >> y4;
        CompG::vec2 v0(x2-x1, y2-y1);
        CompG::vec2 v1(x4-x3, y4-y3);
        double dot = CompG::dot(v0, v1);
        double cross = CompG::cross(v0,v1);
        double absv0 = CompG::abs(v0);
        double absv1 = CompG::abs(v1);
        if (CompG::doubleEqual(dot, absv0*absv1) ||
            CompG::doubleEqual(dot, -absv0*absv1)) {
            // 平行
            std::cout << 2 << "\n";
        } else if (CompG::doubleEqual(dot, 0.0)) {
            // 垂直
            std::cout << 1 << "\n";
        } else {
            std::cout << 0 << "\n";
        }
    }
    return 0;
}