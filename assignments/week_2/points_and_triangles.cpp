// Include header file for standard input/output stream library
#include <iostream>
// Include useful maths operation
#include <cmath>

class Point
{
    public:
        // CONSTRUCTORS
        Point()
        :   x(0),
            y(0)
        {}

        Point(const double &x, const double &y)
        :   x(x),
            y(y)
        {}

        // OPERATORS
        Point operator+(const Point &other) const
        { 
            // create a variable Point to store the result
            Point result (x + other.x, y + other.y);
            return result;
        }

        Point operator-(const Point &other) const 
        { 
            // create a variable Point to store the result
            Point result (x - other.x, y - other.y);
            return result;
        }

        Point &operator+=(const Point &other)
        {
            x = x + other.x;
            y = y + other.y;
            return *this;
        }

        Point &operator-=(const Point &other)
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
        double distance() const
        {
            double distance_origin;
            distance_origin = sqrt(pow(x, 2) + pow(y, 2));
            return distance_origin;
        }

        double distance(const Point &other) const 
        {
            double distance_points;
            distance_points = sqrt(pow((x - other.x), 2) + pow((y - other.y), 2));
            return distance_points;
        }

        Point rotated(double angle) const
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

        Point rotated(double angle, const Point &other) const
        {
            // Initialization of the output
            Point point_rotated;
            // Simplification of the notation
            double s = sin(angle);
            double c = cos(angle);

            // Creating temporary variable that are not constant and can be modified
            double x_temp = x;
            double y_temp = y;
            // Translating the point back to the origin
            x_temp -= other.x;
            y_temp -= other.y;

            // rotate point
            point_rotated.x = x_temp * c - y_temp * s;
            point_rotated.y = x_temp * s + y_temp * c;

            // Applying the translation back again on point_rotated
            point_rotated += other;

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

        Point& rotate(double angle, const Point &other) 
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

            return *this;
        }

};

class Triangle
{
    public:
        // CONSTRUCTORS
        Triangle()
        :   a(),
            b(),
            c()
        {}

        Triangle(const Point &A, const Point &B, const Point &C)
        :   a(A.x, A.y),
            b(B.x, B.y),
            c(C.x, C.y)
        {}

        // OPERATORS    
        // ATTRIBUTES
        Point a, b, c;

        // METHODS
        Triangle translated(const Point &t) const 
        {
            // Create the output variable
            Triangle translated_triang(a, b, c);

            // Applying the translation
            translated_triang.a += t;
            translated_triang.b += t;
            translated_triang.c += t;

            return translated_triang;
        }

        Triangle& translate(const Point &t)
        {
            // Applying the translation
            a += t;
            b += t;
            c += t;

            return *this;
        }

        Triangle rotated(double angle) const
        {
            // creating the new triangle to return
            Triangle rotated_triang = *this;
            
            // rotating each element of the new triangle
            rotated_triang.a.rotate(angle);
            rotated_triang.b.rotate(angle);
            rotated_triang.c.rotate(angle);

            return rotated_triang;
        }

        Triangle rotated(double angle, const Point &other) const
        {
            // creating the new triangle to return
            Triangle rotated_triang = *this;
            
            // rotating each element of the new triangle around the input Point
            rotated_triang.a.rotate(angle, other);
            rotated_triang.b.rotate(angle, other);
            rotated_triang.c.rotate(angle, other);

            return rotated_triang;
        }

        Triangle &rotate(double angle) 
        {
            
            // rotating each element of the triangle
            a.rotate(angle);
            b.rotate(angle);
            c.rotate(angle);

            return *this;
        }

        Triangle &rotate(double angle, const Point &other)
        {
            
            // rotating each element of the triangle around the input Point
            a.rotate(angle, other);
            b.rotate(angle, other);
            c.rotate(angle, other);

            return *this;
        }

        double area() const
        {
            // Using the shoelace formula
            double area = 0.5 * abs(a.x * b.y - c.x * b.y + c.x * a.y - a.x * c.y + b.x * c.y - b.x * a.y);
            return area;
        }
};

int main(){
    // ASSIGNMENT PART I
    Point A(2, 4);
    Point B(3, 5);
    Point C(4, 5);

    B += A;

    double dist_orig = A.distance();
    double dist_points = B.distance(C);

    double pi = M_PI;
    Point D = A.rotated(pi / 2);
    Point E = A.rotated(pi/2, B);

    Point F(0, 1);
    Point G(2, 4);
    F.rotate(pi / 2);
    G.rotate(pi/2, B);

    // Testing Part I of the assignment
    std::cout << "PART I\n" << "A: " << A.x << " " << A.y << std::endl;
    std::cout << "B: " << B.x << " " << B.y << std::endl;
    std::cout << "C: " << C.x << " " << C.y << std::endl;
    std::cout << "distance of A from the origin: " << dist_orig << std::endl;
    std::cout << "distance between B and C: " << dist_points << std::endl;
    std::cout << "Rotation of A around the origin: " << D.x << " " << D.y << std::endl;
    std::cout << "Rotation of A around B: " << E.x << " " << E.y << std::endl;
    std::cout << "Rotation of F around the origin: " << F.x << " " << F.y << std::endl;
    std::cout << "Rotation of G around B: " << E.x << " " << E.y << std::endl;

    // ASSIGNMENT PART II
    Triangle test_1;
    Triangle ABC(A, B, C);
    Triangle ACD(A, C, D);
    Triangle BCD(B, C, D);

    Triangle translated_ABC = ABC.translated(A);
    test_1.translate(C);

    Triangle test_2 = ABC.rotated(pi);
    Triangle test_3 = ABC.rotated(pi, A);

    ACD.rotate(pi);
    BCD.rotate(pi, A);

    // For the area we create a simply triangle to verify the correctness of the code
    Point X; //(0, 0)
    Point Y(0, 3);
    Point Z(5, 0);
    Triangle XYZ(X, Y, Z);
    double area_XYZ = XYZ.area();

    // Testing Part II of the assignment
    std::cout << "\nPART II\n" << "Triangle ABC is: " << "Vertex 1: " << ABC.a.x << " " << ABC.a.y << " Vertex 2: " << ABC.b.x << " " << ABC.b.y << " Vertex 3: " << ABC.c.x << " " << ABC.c.y << std::endl;
    std::cout << "Triangle ABC translated by A is: " << "Vertex 1: " << translated_ABC.a.x << " " << translated_ABC.a.y << " Vertex 2: " << translated_ABC.b.x << " " << translated_ABC.b.y << " Vertex 3: " << translated_ABC.c.x << " " << translated_ABC.c.y << std::endl;
    std::cout << "Triangle test_1 translated by C is: " << "Vertex 1: " << test_1.a.x << " " << test_1.a.y << " Vertex 2: " << test_1.b.x << " " << test_1.b.y << " Vertex 3: " << test_1.c.x << " " << test_1.c.y << std::endl;
    std::cout << "Triangle ABC rotated by pi: " << "Vertex 1: " << test_2.a.x << " " << test_2.a.y << " Vertex 2: " << test_2.b.x << " " << test_2.b.y << " Vertex 3: " << test_2.c.x << " " << test_2.c.y << std::endl;
    std::cout << "Triangle ABC rotated by pi around A: " << "Vertex 1: " << test_3.a.x << " " << test_3.a.y << " Vertex 2: " << test_3.b.x << " " << test_3.b.y << " Vertex 3: " << test_3.c.x << " " << test_3.c.y << std::endl;
    std::cout << "Triangle ACD rotated by pi: " << "Vertex 1: " << ACD.a.x << " " << ACD.a.y << " Vertex 2: " << ACD.b.x << " " << ACD.b.y << " Vertex 3: " << ACD.c.x << " " << ACD.c.y << std::endl;
    std::cout << "Triangle BCD rotated by pi around A: " << "Vertex 1: " << BCD.a.x << " " << BCD.a.y << " Vertex 2: " << BCD.b.x << " " << BCD.b.y << " Vertex 3: " << BCD.c.x << " " << BCD.c.y << std::endl;
    std::cout << "Area of triangle XYZ: " <<  area_XYZ << std::endl;

}