#include <bits/stdc++.h>

namespace Comp2G{
        
    constexpr double EPS = 1e-7;

    struct Line;
    struct Point;
    using vec2 = Point;
    bool doubleEqual(double a, double b); 
    double dot(vec2 &a, vec2 &b);
    struct Point{
        Point(double arg_x,double arg_y)
            :x(arg_x), y(arg_y) {}
        Point(Point a, Point b):x(b.x - a.x), y(b.y - a.y) {}
        Point() {x = y = 0.0;}

        bool operator==(const Point &p) { return (doubleEqual(x, p.x) && doubleEqual(y, p.y));}
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
        Point toUnit() {return Point(x/sqrt(x*x+y*y), y/sqrt(x*x+y*y));}

        using vec2 = Point;
        bool _isHomotropic(vec2 &a, vec2 &b) 
        {
            if (doubleEqual(dot(a, b), a.abs()*b.abs()))
                return true;
            return false;
        }
        bool isHomotropic(vec2 a, vec2 b) { return _isHomotropic(a, b);}
        double x,y;
    };

    using vec2 = Point;
    
    struct Line{
        Line()=default;
        Line(Point a,Point b):p({a,b}){}
        Line(double x1, double y1, double x2, double y2): p({Point(x1,y1),Point(x2,y2)}) {}
        vec2 toVec2() const {return {p[0].x - p[1].x, p[0].y - p[1].y};}
        bool isOnLine(Point a)
        {
            vec2 z = toVec2(); // 1->0
            vec2 x(a, p[0]), y(p[1], a);
            return z == x + y && 
                    a.x <= std::max(p[0].x, p[1].x) &&
                    a.x >= std::min(p[0].x, p[1].x) &&
                    a.y <= std::max(p[0].y, p[1].y) &&
                    a.y >= std::min(p[0].y, p[1].y);
        }
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
    bool doubleEqual(double a, double b) { return fabs(a - b) <= EPS; }

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
        double dot = Comp2G::dot(a, b), absa = abs(a), absb = abs(b);
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
        double dot = Comp2G::dot(a, b);
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
        return crossVal < 0 ? true : false;
    }
    bool transverse(Line a, Line b) {return _transverse(a, b);}

    // 如果两线段相交则返回true
    bool _isIntersection(Line &a, Line &b)
    {
        if (_quickRejection(a, b))
            if (_transverse(a, b) && _transverse(b, a))
                return true;
        return false;
    }
    bool isIntersection(Line a, Line b) {return _isIntersection(a, b);}

    // 返回两条线段/或直线的交点，使用前需先判断直线相交 
    Point _CrossPoint(Line &a, Line &b)
    {
        vec2 vu = b.p[0] - a.p[0]; // b0->a0
        vec2 va = a.toVec2(), vb = b.toVec2();
        double i = std::abs(cross(va, vb)) * std::abs(cross(vu, vb));
        vec2 Ta = va.toUnit()*i;
        return a.p[0] + Ta;
    }

    double _distance(Point &a, Point &b) { double x = a.x-b.x, y = a.y-b.y; return sqrt(x*x+y*y);}
    double distance(Point a, Point b) {return _distance(a, b);}
};

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    int t;
    std::cin >> t;
    while (t--) {
        double x1,x2,x3,x4,y1,y2,y3,y4;
        std::cin >> x1 >> y1;
        std::cin >> x2 >> y2;
        std::cin >> x3 >> y3;
        std::cin >> x4 >> y4;
        Comp2G::Line a(x1, y1, x2, y2), b(x3, y3, x4, y4);
        double ans = .0;
        if (Comp2G::_isIntersection(a, b)) ;
        else {
            ans = 1e18;
            Comp2G::Point pa = Comp2G::projection(a, b.p[0]);
            if (a.isOnLine(pa))
                ans = std::min(ans, Comp2G::distance(pa, b.p[0]));
            pa = Comp2G::projection(a, b.p[1]);
            if (a.isOnLine(pa))
                ans = std::min(ans, Comp2G::distance(pa, b.p[1]));
            pa = Comp2G::projection(b, a.p[0]);
            if (b.isOnLine(pa))
                ans = std::min(ans, Comp2G::distance(pa, a.p[0]));
            pa = Comp2G::projection(b, a.p[1]);
            if (b.isOnLine(pa))
                ans = std::min(ans, Comp2G::distance(pa, a.p[1]));
            ans = std::min(ans, Comp2G::distance(a.p[0], b.p[0]));
            ans = std::min(ans, Comp2G::distance(a.p[0], b.p[1]));
            ans = std::min(ans, Comp2G::distance(a.p[1], b.p[0]));
            ans = std::min(ans, Comp2G::distance(a.p[1], b.p[1]));
        }
        std::cout << std::setprecision(10) << std::fixed << ans << "\n";
    }
    return 0;

}