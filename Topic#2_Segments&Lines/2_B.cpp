
#include <bits/stdc++.h>

namespace CompG{
        
    constexpr double EPS = 1e-7;

    struct Line;

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

        double x,y;
    };

    using vec2 = Point;
    
    struct Line{
        Line()=default;
        Line(Point a,Point b):p({a,b}){}
        Line(double x1, double y1, double x2, double y2): p({Point(x1,y1),Point(x2,y2)}) {}
        vec2 toVec2() const {return {p[0].x - p[1].x, p[0].y - p[1].y};}
        /*void outPut(){
            printf("(%.3lf,%.3lf)-(%.3lf,%.3lf)\n",p[0].x,p[0].y,p[1].x,p[1].y);
        }*/
        Point p[2];
    };

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

    // 判断两向量是否为平行向量
    bool _isParallel(vec2 &a, vec2 &b)
    {
        double dot = CompG::dot(a, b), absa = abs(a), absb = abs(b);
        return (doubleEqual(dot, absa*absb) || doubleEqual(dot, -absa*absb)) ? 1 : 0;
    }

    // 判断两向量是否为平行向量
    bool isParallel(vec2 a, vec2 b)
    {
        return _isParallel(a, b);
    }

    // 两向量垂直
    bool _isOrthogonal(vec2 &a, vec2 &b)
    {
        double dot = CompG::dot(a, b);
        return doubleEqual(dot, .0) ? 1 : 0;
    } 

    // 两向量垂直
    bool isOrthogonal(vec2 a, vec2 b) { return _isOrthogonal(a, b); }

    bool _QJcompare(double x1, double x2, double x3, double x4) 
    {
        if (x1 < x3 && x2 < x3 && x1 < x4 && x2 < x4) return true;
        if (x3 < x1 && x4 < x1 && x3 < x2 && x4 < x2) return true;
        return false;
    }

    bool _quickRejection(Line &a, Line &b)
    {
        if (_QJcompare(a.p[0].x, a.p[1].x, b.p[0].x, b.p[1].x) ||
            _QJcompare(a.p[0].y, a.p[1].y, b.p[0].y, b.p[1].y)) return false;
        else return true;
    }
    bool quickRejection(Line a, Line b) {return _quickRejection(a, b);}

    bool _transverse(Line &a, Line &b) 
    {
        vec2 va = a.toVec2();
        vec2 ab1(b.p[0] - a.p[0]), ab2(b.p[1] - a.p[0]);
        double crossVal = cross(va, ab1) * cross(va, ab2);
        return crossVal <= 0 ? true : false;
    }
    bool transverse(Line a, Line b) {return _transverse(a, b);}
};


signed main()
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
        CompG::Line a(x1, y1, x2, y2), b(x3, y3, x4, y4);
        int ans = 0;
        // if (CompG::isParallel(a.toVec2(), b.toVec2())) 
        //     ans = 0;
        // else {
        //     if (CompG::_quickRejection(a, b))
        //         if (CompG::_transverse(a, b) && CompG::_transverse(b, a))
        //             ans = 1;
        // }
        if (CompG::_quickRejection(a, b))
            if (CompG::_transverse(a, b) && CompG::_transverse(b, a))
                ans = 1;
        std::cout << ans << "\n";
    }    
    return 0;
}