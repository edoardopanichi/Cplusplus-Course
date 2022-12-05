// Include header file for standard input/output stream library
#include <iostream>
// Include useful maths operation
#include <cmath>
#include <numbers>

class Point
{
    public:
        // CONSTRUCTORS
        Point()
        :   x(0),
            y(0)
        {}

        Point(double x, double y)
        :   x(x),
            y(y)
        {}

        // OPERATORS
        Point operator+(Point other) 
        { 
            // create a variable Point to store the result
            Point result (x + other.x, y + other.y);
            return result;
        }

        Point operator-(Point other) 
        { 
            // create a variable Point to store the result
            Point result (x - other.x, y - other.y);
            return result;
        }

        Point &operator+=(Point other) 
        {
            x = x + other.x;
            y = y + other.y;
            return *this;
        }

        Point &operator-=(Point other) 
        { 
            x = x - other.x;
            y = y - other.y;
            return *this;
        }

        // ATTRIBUTES
        double x;
        double y;

        // METHODS
        // this function returns the distance from the origin
        double distance() 
        {
            double distance_origin;
            distance_origin = sqrt(pow(x, 2) + pow(y, 2));
            return distance_origin;
        }

        double distance(Point other) 
        {
            double distance_points;
            distance_points = sqrt(pow((x - other.x), 2) + pow((y - other.y), 2));
            return distance_points;
        }

        Point rotated(double angle)
        {
            // Initialization of the output
            Point point_rotated;
            // Simplification of the notation
            double s = sin(angle);
            double c = cos(angle);

            // rotate point
            point_rotated.x = x * c - y * s;
            point_rotated.y = x * s + y * c;

            return point_rotated;
        }

        Point rotated(double angle, Point other)
        {
            // Initialization of the output
            Point point_rotated;
            // Simplification of the notation
            double s = sin(angle);
            double c = cos(angle);

            // Translating the point back to the origin
            x -= other.x;
            y -= other.y;

            // rotate point
            point_rotated.x = x * c - y * s;
            point_rotated.y = x * s + y * c;

            // Applying the translation back again (this time also on point_rotated)
            point_rotated += other;
            x += other.x;
            y += other.y;

            return point_rotated;
        }

        Point& rotate(double angle) 
        { 
            // Simplification of the notation
            double s = sin(angle);
            double c = cos(angle);

            double x_temp = x;

            // rotate point
            x = x * c - y * s;
            y = x_temp * s + y * c;

            return *this; 
        }

        Point& rotate(double angle, const Point other) 
        {
            // Simplification of the notation
            double s = sin(angle);
            double c = cos(angle);

            // Translating the point back to the origin
            x -= other.x;
            y -= other.y;

            double x_temp = x;

            // rotate point
            x = x * c - y * s;
            y = x_temp * s + y * c;

            // Applying the translation back again 
            x += other.x;
            y += other.y;
            this -> x += other.y;

            return *this;
        }

};

int main(){
    Point A(2, 4);
    Point B = {3, 5};
    Point C(4, 5);

    B += A;

    double dist_orig = A.distance();
    double dist_points = B.distance(C);

    double pi = atan(1)*4;
    Point D = A.rotated(pi/2);
    Point E = A.rotated(pi/2, B);

    Point F = {0, 1};
    Point G = {2, 4};
    F.rotate(pi / 2);
    G.rotate(pi/2, B);


    // Testing Part I of the assignment
    std::cout << "A: " << A.x << " " << A.y << std::endl;
    std::cout << "B: " << B.x << " " << B.y << std::endl;
    std::cout << "C: " << C.x << " " << C.y << std::endl;
    std::cout << "distance of A from the origin: " << dist_orig << std::endl;
    std::cout << "distance between B and C: " << dist_points << std::endl;
    std::cout << "Rotation of A around the origin: " << D.x << " " << D.y << std::endl;
    std::cout << "Rotation of A around B: " << E.x << " " << E.y << std::endl;
    std::cout << "Rotation of F around the origin: " << F.x << " " << F.y << std::endl;
    std::cout << "Rotation of G around B: " << E.x << " " << E.y << std::endl;

}