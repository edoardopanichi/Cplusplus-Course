#include <functional>
#include <iostream>
#include <cmath>

// Your implementation of the abstract base class goes here
class Derivative
{
public: //[DO NOT MODIFY/REMOVE THIS GETTER: IT IS USED IN THE SPECTEST]
    double GetH() const {return h;}

    // Constructors
    // Default constructor
    Derivative() 
    : h(1e-8) {}
    // Constructor with h as parameter
    Derivative(double h)
    : h(h) {}

    // Implementation of the pure virtual function goes here
    virtual double differentiate(const std::function<double(double)>& func, double x) const = 0;

protected:
    double h;
};

// Your implementation of the derived class for the central difference scheme goes here
class CentralDifference: public Derivative
{
public:
    // Constructors (use the base class constructors)
    using Derivative::Derivative;

    // Implementation of the pure virtual function goes here
    double differentiate(const std::function<double(double)>& func, double x) const override
    {
        return (func(x+h) - func(x-h)) / (2*h);
    }
};

// Your implementation of the derived class for the forward difference scheme goes here
class ForwardDifference : public Derivative
{
    public:
    // Constructors (use the base class constructors)
    using Derivative::Derivative;

    // Implementation of the pure virtual function goes here
    double differentiate(const std::function<double(double)>& func, double x) const override
    {
        return (func(x+h) - func(x)) / h;
    }
};

// The implementation of the traditional function goes here
const double myfunc_2(double x) 
{
    return (x * x * x);
}

int main()
{
    // Your tests go here

    // First two tests are with lambda functions
    CentralDifference test;
    double derivative_test = test.differentiate([](double x){return x*x;}, 3.0);
    std::cout << derivative_test << std::endl;

    ForwardDifference test2;
    double derivative_test2 = test2.differentiate([](double x){return x*x;}, 3.0);
    std::cout << derivative_test2 << std::endl;

    // Last two tests are with traditional functions
    CentralDifference test3;
    double derivative_test3 = test3.differentiate(myfunc_2, 3.0);
    std::cout << derivative_test3 << std::endl;

    ForwardDifference test4;
    double derivative_test4 = test4.differentiate(myfunc_2, 3.0);
    std::cout << derivative_test4 << std::endl;
    

    return 0;
}