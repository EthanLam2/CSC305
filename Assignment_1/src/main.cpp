////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    // TODO
    return u.x() * v.y() - u.y() * v.x();
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    // TODO

    double det1 = det(b-a,c-a);
    double det2 = det(b-a,d-a);
    double det3 = det(d-c,a-c);
    double det4 = det(d-c,b-c);
    double prod1 = det1*det2;
    double prod2 = det3*det4;

    return (prod1 < 0 && prod2 < 0);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    // TODO
    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    // TODO
    Array2d minValues;
    Array2d maxValues;
    minValues.setConstant(std::numeric_limits<double>::infinity());
    maxValues.setConstant(-std::numeric_limits<double>::infinity());

    for (const Vector2d &point : poly) {
        minValues = minValues.min(point.array());
        maxValues = maxValues.max(point.array());
    }
    Vector2d outside(maxValues[0] + 1, query.y()); 
    int intersectionCount = 0;

    for (size_t i = 0; i < poly.size(); ++i) {
        const Vector2d &p1 = poly[i];
        const Vector2d &p2 = poly[(i + 1) % poly.size()];

       
        if (intersect_segment(p1, p2, query, outside)) {
            ++intersectionCount;
        }
    }

    return (intersectionCount % 2 == 1);
    
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);
    // TODO
    int num_points;
    in >> num_points;
    double x;
    double y;
    double z;
    for (int i=0; i < num_points; i++)
    {
        in >> x >> y >> z;
        points.push_back(Vector2d(x,y));
    }
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out (filename);
    // TODO
    out << points.size() << '\n';

    for (const auto &point : points)
    {
        out << point.x() << ' ' << point.y() << " 0\n";
    }

}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    
    
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);
    
    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}