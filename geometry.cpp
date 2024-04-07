#include <bits/stdc++.h>
using namespace std;
typedef double T;
struct pt
{
    T x, y;
    pt operator+(pt p) { return {x + p.x, y + p.y}; }
    pt operator-(pt p) { return {x - p.x, y - p.y}; }
    pt operator*(T d) { return {x * d, y * d}; }
    pt operator/(T d) { return {x / d, y / d}; } // only for floatingpoint
};
bool operator==(pt a, pt b) { return a.x == b.x && a.y == b.y; }
bool operator!=(pt a, pt b) { return !(a == b); }
T sq(pt p) { return p.x * p.x + p.y * p.y; }
double abs(pt p) { return sqrt(sq(p)); }
int sgn(T x) { return (T(0) < x) - (x < T(0)); }
pt translate(pt v, pt p) { return p + v; }
pt scale(pt c, double factor, pt p) { return c + (p - c) * factor; }
pt rot(pt p, double a) { return {p.x * cos(a) - p.y * sin(a), p.x * sin(a) + p.y * cos(a)}; }
pt perp(pt p) { return {-p.y, p.x}; }
T dot(pt v, pt w) { return v.x * w.x + v.y * w.y; }
double angle(pt v, pt w)
{
    double cosTheta = dot(v, w) / abs(v) / abs(w);
    return acos(max(-1.0, min(1.0, cosTheta)));
}
T cross(pt v, pt w) { return v.x * w.y - v.y * w.x; }
T orient(pt a, pt b, pt c) { return cross(b - a, c - a); }
bool inAngle(pt a, pt b, pt c, pt p)
{
    assert(orient(a, b, c) != 0);
    if (orient(a, b, c) < 0)
        swap(b, c);
    return orient(a, b, p) >= 0 && orient(a, c, p) <= 0;
}
bool isConvex(vector<pt> p)
{
    bool hasPos = false, hasNeg = false;
    for (int i = 0, n = p.size(); i < n; i++)
    {
        int o = orient(p[i], p[(i + 1) % n], p[(i + 2) % n]);
        if (o > 0)
            hasPos = true;
        if (o < 0)
            hasNeg = true;
    }
    return !(hasPos && hasNeg);
}

// line
struct line
{
    pt v;
    T c;
    // From direction vector v and offset c
    line(pt v, T c) : v(v), c(c) {}
    // From equation ax+by=c
    line(T a, T b, T c) : v({b, -a}), c(c) {}
    // From points P and Q
    line(pt p, pt q) : v(q - p), c(cross(v, p)) {}
    // Will be defined later:
    // - these work with T = int
    T side(pt p);
    double dist(pt p);
    line perpThrough(pt p);
    bool cmpProj(pt p, pt q);
    line translate(pt t);
    // - these require T = double
    line shiftLeft(double dist);
    pt proj(pt p);
    pt refl(pt p);
};
T line::side(pt p) { return cross(v, p) - c; }
double line::dist(pt p) { return abs(side(p)) / abs(v); }
line line::perpThrough(pt p) { return {p, p + perp(v)}; }
bool line::cmpProj(pt p, pt q) { return dot(v, p) < dot(v, q); }
line line::translate(pt t) { return {v, c + cross(v, t)}; }
line line::shiftLeft(double dist) { return {v, c + dist * abs(v)}; }
bool intersection(line l1, line l2, pt &out)
{
    T d = cross(l1.v, l2.v);
    if (d == 0)
        return false;
    out = (l2.v * l1.c - l1.v * l2.c) / d; // requires floating-pointcoordinates
    return true;
}
pt line::proj(pt p) { return p - perp(v) * side(p) / sq(v); }
pt line::refl(pt p) { return p - perp(v) * 2 * side(p) / sq(v); }
line bisector(line l1, line l2, bool interior)
{
    assert(cross(l1.v, l2.v) != 0); // l1 and l2 cannot be parallel!
    double sign = interior ? 1 : -1;
    return {l2.v / abs(l2.v) + l1.v / abs(l1.v) * sign,
            l2.c / abs(l2.v) + l1.c / abs(l1.v) * sign};
}
bool inDisk(pt a, pt b, pt p) { return dot(a - p, b - p) <= 0; }//if p is indide circle of diameter ab
bool onSegment(pt a, pt b, pt p) { return orient(a, b, p) == 0 && inDisk(a, b, p); }

// segment segment intersection
// case: 01 (interior)
bool properInter(pt a, pt b, pt c, pt d, pt &out)
{
    double oa = orient(c, d, a),
           ob = orient(c, d, b),
           oc = orient(a, b, c),
           od = orient(a, b, d);
    // Proper intersection exists iff opposite signs
    if (oa * ob < 0 && oc * od < 0)
    {
        out = (a * ob - b * oa) / (ob - oa);
        return true;
    }
    return false;
}
// case: 02 (corners)
//  To create sets of points we need a comparison function
struct cmpX
{
    bool operator()(pt a, pt b)
    {
        return make_pair(a.x, a.y) < make_pair(b.x, b.y);
    }
};
set<pt, cmpX> inters(pt a, pt b, pt c, pt d)
{
    pt out;
    if (properInter(a, b, c, d, out))
        return {out};
    set<pt, cmpX> s;
    if (onSegment(c, d, a))
        s.insert(a);
    if (onSegment(c, d, b))
        s.insert(b);
    if (onSegment(a, b, c))
        s.insert(c);
    if (onSegment(a, b, d))
        s.insert(d);
    return s;
}

// Segment-point distance
double segPoint(pt a, pt b, pt p)
{
    if (a != b)
    {
        line l(a, b);
        if (l.cmpProj(a, p) && l.cmpProj(p, b)) // if closest to projection
            return l.dist(p);                   // output distance toline
    }
    return min(abs(p - a), abs(p - b)); // otherwise distance to A or B
}

// segment-segment distance
double segSeg(pt a, pt b, pt c, pt d)
{
    pt dummy;
    if (properInter(a, b, c, d, dummy))
        return 0;
    return min({segPoint(a, b, c), segPoint(a, b, d),
                segPoint(c, d, a), segPoint(c, d, b)});
}

// area of triangle & polygon
double areaTriangle(pt a, pt b, pt c) { return abs(cross(b - a, c - a)) / 2.0; }
double areaPolygon(vector<pt> p)
{
    double area = 0.0;
    for (int i = 0, n = p.size(); i < n; i++)
    {
        area += cross(p[i], p[(i + 1) % n]); // wrap back to 0 if i == n-1
    }
    return abs(area) / 2.0;
}

// Point inside polygon (Cutting-ray test)

//  true if P at least as high as A (blue part)
bool above(pt a, pt p) { return p.y >= a.y; }
// check if [PQ] crosses ray from A
bool crossesRay(pt a, pt p, pt q) { return (above(a, q) - above(a, p)) * orient(a, p, q) > 0; }
// if strict, returns false when A is on the boundary
bool inPolygon(vector<pt> p, pt a, bool strict = true)
{
    int numCrossings = 0;
    for (int i = 0, n = p.size(); i < n; i++)
    {
        if (onSegment(p[i], p[(i + 1) % n], a))
            return !strict;
        numCrossings += crossesRay(a, p[i], p[(i + 1) % n]);
    }
    return numCrossings & 1; // inside if odd number of crossings
}

// winding number
// amplitude travelled around point A, from P to Q
bool half(pt p) { return p.y > 0 || (p.y == 0 && p.x < 0); }
struct angle
{
    pt d;
    int t = 0;    // direction and number of full turns
    angle t180() { return {d * (-1), t + half(d)}; }
    angle t360() { return {d, t + 1}; }
    bool operator<(const angle &b) const
    {
        return make_tuple(t, half(d), 0) <
               make_tuple(b.t, half(b.d), cross(d, b.d));
    }
};

angle moveTo(angle a, pt newD)
{
    // check that segment [DD’] doesn’t go through the origin
    assert(!onSegment(a.d, newD, {0, 0}));
    angle b{newD, a.t};
    if (a.t180() < b) // if b more than half a turn bigger
        b.t--;        // decrease b by a full turn
    if (b.t180() < a) // if b more than half a turn smaller
        b.t++;        // increase b by a full turn
    return b;
}
int windingNumber(vector<pt> p, pt a)
{
    angle a{p.back()}; // start at last vertex
    for (pt d : p)
        a = moveTo(a, d); // move to first vertex, second, etc.
    return a.t;
}

// circle

// center
pt circumCenter(pt a, pt b, pt c)
{
    b = b - a, c = c - a;     // consider coordinates relative to A
    assert(cross(b, c) != 0); // no circumcircle if A,B,C aligned
    return a + perp(b * sq(c) - c * sq(b)) / cross(b, c) / 2;
}

// circle-line intersection
int circleLine(pt o, double r, line l, pair<pt, pt> &out)
{
    double h2 = r * r - l.sqDist(o);
    if (h2 >= 0)
    {                                     // the line touches the circle
        pt p = l.proj(o);                 // point P
        pt h = l.v * sqrt(h2) / abs(l.v); // vector parallel to l, of
        length h
            out = {p - h, p + h};
    }
    return 1 + sgn(h2);
}

// circle-circle intersection
int circleCircle(pt o1, double r1, pt o2, double r2, pair<pt, pt> &out)
{
    pt d = o2 - o1;
    double d2 = sq(d);
    if (d2 == 0)
    {
        assert(r1 != r2);
        return 0;
    }                                         // concentric circles
    double pd = (d2 + r1 * r1 - r2 * r2) / 2; // = |O_1P| * d
    double h2 = r1 * r1 - pd * pd / d2;       // = hˆ2
    if (h2 >= 0)
    {
        pt p = o1 + d * pd / d2, h = perp(d) * sqrt(h2 / d2);
        out = {p - h, p + h};
    }
    return 1 + sgn(h2);
}

// Tangents of circle

/*
the same code can be used to find the tangent to a circle
passing through a point by setting r2 to 0 (in which case the value of inner
doesn’t matter).
*/
int tangents(pt o1, double r1, pt o2, double r2, bool inner, vector<pair<pt, pt>> &out)
{
    if (inner)
        r2 = -r2;
    pt d = o2 - o1;
    double dr = r1 - r2, d2 = sq(d), h2 = d2 - dr * dr;
    if (d2 == 0 || h2 < 0)
    {
        assert(h2 != 0);
        return 0;
    }
    for (double sign : {-1, 1})
    {
        pt v = (d * dr + perp(d) * sqrt(h2) * sign) / d2;
        out.push_back({o1 + v * r1, o2 + v * r2});
    }
    return 1 + (h2 > 0);
}

// 3D Geometry

typedef double T;
struct p3
{
    T x, y, z;
    // Basic vector operations
    p3 operator+(p3 p) { return {x + p.x, y + p.y, z + p.z}; }
    p3 operator-(p3 p) { return {x - p.x, y - p.y, z - p.z}; }
    p3 operator*(T d) { return {x * d, y * d, z * d}; }
    p3 operator/(T d) { return {x / d, y / d, z / d}; } // only for floating - point
    // Some comparators
    bool
    operator==(p3 p)
    {
        return tie(x, y, z) == tie(p.x, p.y, p.z);
    }
    bool operator!=(p3 p) { return !operator==(p); }
    bool operator<(const p3 &q) const { return tie(x, y, z) < tie(q.x, q.y, q.z); }
};
// smaller angle between two vector
double smallAngle(p3 v, p3 w) { return acos(min(abs(v | w) / abs(v) / abs(w), 1.0)); }
// dot product
T operator|(p3 v, p3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }
T sq(p3 v) { return v | v; }
double abs(p3 v) { return sqrt(sq(v)); }
p3 unit(p3 v) { return v / abs(v); }
double angle(p3 v, p3 w)
{
    double cosTheta = (v | w) / abs(v) / abs(w);
    return acos(max(-1.0, min(1.0, cosTheta)));
}

// cross product
p3 operator*(p3 v, p3 w)
{
    return {v.y * w.z - v.z * w.y,
            v.z * w.x - v.x * w.z,
            v.x * w.y - v.y * w.x};
}
/*
It is positive if S is on the side of plane P QR in the direction of PQ×PR,
negative if S is on the other side, and zero if S is on the plane
*/
T orient(p3 p, p3 q, p3 r, p3 s) { return (q - p) * (r - p) | (s - p); }
T orientByNormal(p3 p, p3 q, p3 r, p3 n) { return (q - p) * (r - p) | n; }

struct plane
{
    p3 n;
    T d;
    // From normal n and offset d
    plane(p3 n, T d) : n(n), d(d) {}
    // From normal n and point P
    plane(p3 n, p3 p) : n(n), d(n | p) {}
    // From three non-collinear points P,Q,R
    plane(p3 p, p3 q, p3 r) : plane((q - p) * (r - p), p) {}
    // Will be defined later:
    // - these work with T = int
    T side(p3 p);
    double dist(p3 p);
    plane translate(p3 t);
    // - these require T = double
    plane shiftUp(double dist);
    p3 proj(p3 p);
    p3 refl(p3 p);
};

T plane::side(p3 p) { return (n | p) - d; }
double plane::dist(p3 p) { return abs(side(p)) / abs(n); }
plane plane::translate(p3 t) { return {n, d + (n | t)}; }
plane plane::shiftUp(double dist) { return {n, d + dist * abs(n)}; }
p3 plane::proj(p3 p) { return p - n * side(p) / sq(n); }
p3 plane::refl(p3 p) { return p - n * 2 * side(p) / sq(n); }

// 3d to 2d conoversion
struct coords
{
    p3 o, dx, dy, dz;
    // From three points P,Q,R on the plane:
    // build an orthonormal 3D basis
    coords(p3 p, p3 q, p3 r) : o(p)
    {
        dx = unit(q - p);
        dz = unit(dx * (r - p));
        dy = dz * dx;
    }
    // From four points P,Q,R,S:
    // take directions PQ, PR, PS as is
    coords(p3 p, p3 q, p3 r, p3 s) : o(p), dx(q - p), dy(r - p), dz(s - p) {}
    pt pos2d(p3 p) { return {(p - o) | dx, (p - o) | dy}; }
    p3 pos3d(p3 p) { return {(p - o) | dx, (p - o) | dy, (p - o) | dz}; }
};

struct line3d
{
    p3 d, o;
    // From two points P, Q
    line3d(p3 p, p3 q) : d(q - p), o(p) {}
    // From two planes p1, p2 (requires T = double)
    line3d(plane p1, plane p2); // will be implemented later
    // Will be defined later:
    // - these work with T = int
    double sqDist(p3 p);
    double dist(p3 p);
    bool cmpProj(p3 p, p3 q);
    // - these require T = double
    p3 proj(p3 p);
    p3 refl(p3 p);
    p3 inter(plane p);
};
double line3d::sqDist(p3 p) { return sq(d * (p - o)) / sq(d); }
double line3d::dist(p3 p) { return sqrt(sqDist(p)); }
bool line3d::cmpProj(p3 p, p3 q) { return (d | p) < (d | q); }
p3 line3d::proj(p3 p) { return o + d * (d | (p - o)) / sq(d); }
p3 line3d::refl(p3 p) { return proj(p) * 2 - p; }

// Plane line intersection
p3 line3d::inter(plane p) { return o - d * p.side(o) / (d | p.n); }
// Note that this is undefined when n·d = 0

// Plane Plane intersection
line3d(plane p1, plane p2)
{
    d = p1.n * p2.n; // undefined if d=0
    o = (p2.n * p1.d - p1.n * p2.d) * d / sq(d);
}

// line line distance nearest distance
double line3d::dist(line3d l1, line3d l2)
{
    p3 n = l1.d * l2.d;
    if (n == zero) // parallel
        return l1.dist(l2.o);
    return abs((l2.o - l1.o) | n) / abs(n);
}

// nine line closest point
// If we want C2 instead we just have to swap l1 and l2.
p3 closestOnL1(line3d l1, line3d l2)
{
    p3 n2 = l2.d * (l1.d * l2.d);
    return l1.o + l1.d * ((l2.o - l1.o) | n2) / (l1.d | n2);
}

// Algle between two plane
double angle(plane p1, plane p2)
{
    return smallAngle(p1.n, p2.n);
}
bool isParallel(plane p1, plane p2)
{
    return p1.n * p2.n == zero;
}
bool isPerpendicular(plane p1, plane p2)
{
    return (p1.n | p2.n) == 0;
}

// Angle between two line
double angle(line3d l1, line3d l2)
{
    return smallAngle(l1.d, l2.d);
}
bool isParallel(line3d l1, line3d l2)
{
    return l1.d * l2.d == zero;
}
bool isPerpendicular(line3d l1, line3d l2)
{
    return (l1.d | l2.d) == 0;
}

// Angle between a plane and a line
double angle(plane p, line3d l)
{
    return M_PI / 2 - smallAngle(p.n, l.d);
}
bool isParallel(plane p, line3d l)
{
    return (p.n | l.d) == 0;
}
bool isPerpendicular(plane p, line3d l)
{
    return p.n * l.d == zero;
}

line3d perpThrough(plane p, p3 o) { return line(o, o + p.n); }
plane perpThrough(line3d l, p3 o) { return plane(l.d, o); }

// surface area of a polygon in 3D
p3 vectorArea2(vector<p3> p)
{ // vector area * 2 (to avoiddivisions)
    p3 S = zero;
    for (int i = 0, n = p.size(); i < n; i++)
        S = S + p[i] * p[(i + 1) % n];
    return S;
}
double area(vector<p3> p) { return abs(vectorArea2(p)) / 2.0; }

struct edge
{
    int v;
    bool same; // = is the common edge in the same order?
};
// Given a series of faces (lists of points), reverse some of them
// so that their orientations are consistent
void reorient(vector<vector<p3>> &fs)
{
    int n = fs.size();
    // Find the common edges and create the resulting graph
    vector<vector<edge>> g(n);
    map<pair<p3, p3>, int> es;
    for (int u = 0; u < n; u++)
    {
        for (int i = 0, m = fs[u].size(); i < m; i++)
        {
            p3 a = fs[u][i], b = fs[u][(i + 1) % m];
            // Let’s look at edge [AB]
            if (es.count({a, b}))
            { // seen in same order
                int v = es[{a, b}];
                g[u].push_back({v, true});
                g[v].push_back({u, true});
            }
            else if (es.count({b, a}))
            { // seen in different order
                int v = es[{b, a}];
                g[u].push_back({v, false});
                g[v].push_back({u, false});
            }
            else
            { // not seen yet
                es[{a, b}] = u;
            }
        }
    }
    // Perform BFS to find which faces should be flipped
    vector<bool> vis(n, false), flip(n);
    flip[0] = false;
    queue<int> q;
    q.push(0);
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (auto e : g[u])
        {
            if (!vis[e.v])
            {
                vis[e.v] = true;
                // If the edge was in the same order,
                // exactly one of the two should be flipped
                flip[e.v] = (flip[u] ^ e.same);
                q.push(e.v);
            }
        }
    }
    // Actually perform the flips
    for (int u = 0; u < n; u++)
        if (flip[u])
            reverse(fs[u].begin(), fs[u].end());
}

// volume of polyhedron
double volume(vector<vector<p3>> fs)
{
    double vol6 = 0.0;
    for (vector<p3> f : fs)
        vol6 += (vectorArea2(f) | f[0]);
    return abs(vol6) / 6.0;
}

// Sphere
p3 sph(double r, double lat, double lon) // lat = okkhorekha
{
    lat *= M_PI / 180, lon *= M_PI / 180;
    return {r * cos(lat) * cos(lon), r * cos(lat) * sin(lon), r * sin(lat)};
}

// Shpere line intersection
int sphereLine(p3 o, double r, line3d l, pair<p3, p3> &out)
{
    double h2 = r * r - l.sqDist(o);
    if (h2 < 0)
        return 0;                     // the line doesn’t touch the sphere
    p3 p = l.proj(o);                 // point P
    p3 h = l.d * sqrt(h2) / abs(l.d); // vector parallel to l, of length
    h
        out = {p - h, p + h};
    return 1 + (h2 > 0);
}

// great-circle distance
double greatCircleDist(p3 o, double r, p3 a, p3 b) { return r * angle(a - o, b - o); }

// Spherical segment intersection
bool validSegment(p3 a, p3 b) { return a * b != zero || (a | b) > 0; }
bool properInter(p3 a, p3 b, p3 c, p3 d, p3 &out)
{
    p3 ab = a * b, cd = c * d; // normals of planes OAB and OCD
    int oa = sgn(cd | a),
        ob = sgn(cd | b),
        oc = sgn(ab | c),
        od = sgn(ab | d);
    out = ab * cd * od; // four multiplications => careful with overflow!
    return (oa != ob && oc != od && oa != oc);
}
bool onSphSegment(p3 a, p3 b, p3 p)
{
    p3 n = a * b;
    if (n == zero)
        return a * p == zero && (a | p) > 0;
    return (n | p) == 0 && (n | a * p) >= 0 && (n | b * p) <= 0;
}
struct directionSet : vector<p3>
{
    using vector::vector; // import constructors
    void insert(p3 p)
    {
        for (p3 q : *this)
            if (p * q == zero)
                return;
        push_back(p);
    }
};
directionSet intersSph(p3 a, p3 b, p3 c, p3 d)
{
    assert(validSegment(a, b) && validSegment(c, d));
    p3 out;
    if (properInter(a, b, c, d, out))
        return {out};
    directionSet s;
    if (onSphSegment(c, d, a))
        s.insert(a);
    if (onSphSegment(c, d, b))
        s.insert(b);
    if (onSphSegment(a, b, c))
        s.insert(c);
    if (onSphSegment(a, b, d))
        s.insert(d);
    return s;
}

// Angles on a sphere
double angleSph(p3 a, p3 b, p3 c) { return angle(a * b, a * c); }
double orientedAngleSph(p3 a, p3 b, p3 c)
{
    if ((a * b | c) >= 0)
        return angleSph(a, b, c);
    else
        return 2 * M_PI - angleSph(a, b, c);
}

// Spherical polygons and area
double areaOnSphere(double r, vector<p3> p)
{
    int n = p.size();
    double sum = -(n - 2) * M_PI;
    for (int i = 0; i < n; i++)
        sum += orientedAngleSph(p[(i + 1) % n], p[(i + 2) % n], p[i]);
    return r * r * sum;
}

// solid angle
// if we find the area of a polygon on a sphere making r = 1
// or if we devide the total area by r*r

// winding Number 3D
int windingNumber3D(vector<vector<p3>> fs)
{
    double sum = 0;
    for (vector<p3> f : fs)
        sum += remainder(areaOnSphere(1, f), 4 * M_PI);
    return round(sum / (4 * M_PI));
}
