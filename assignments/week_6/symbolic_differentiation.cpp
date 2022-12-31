#include <iostream>
#include <cmath>

// Implementing ConstInt struct with eval() and derivative methods
template <int value>
struct ConstInt 
{
    // Since the struct holds constant values, eval() returns the value itself
    static constexpr int eval(double x) { return value; }
    // Derivative of a constant value is zero
    using derivative = ConstInt<0>;
};

// This struct takes two struct as template and adds them returning an int value
template <typename struct1, typename struct2>
struct Add
{
    static int eval(double x) { return struct1::eval(x) + struct2::eval(x); }
    typedef Add<typename struct1::derivative, typename struct2::derivative> derivative;
};

// This struct takes two struct as template and multiply them
template <typename struct1, typename struct2>
struct Mul
{
    // Multiplying the values of the two structs
    static constexpr int eval(double x) { return struct1::eval(x) * struct2::eval(x); }
    typedef Add<Mul<struct1, typename struct2::derivative>, Mul<typename struct1::derivative, struct2>> derivative;
};

// This struct represents a monomial function
template <int exponent>
struct Monomial 
{
  static double eval(double x) { return std::pow(x, exponent); }

  typedef Mul<ConstInt<exponent>, Monomial<exponent - 1>> derivative;
};

// Template specialization for monomial function with exponent 0
template <>
struct Monomial<0> 
{
  static double eval(double x) { return 1.0; }
  typedef ConstInt<0> derivative;
};

// Neg struct, eval() returns the negative value of the struct for x
template <typename struct1>
struct Neg
{
    static double eval(double x) { return -struct1::eval(x); }
    typedef Neg<typename struct1::derivative> derivative;
};

// Declaring structs for sine and cosine functions
struct Cos;
struct Sin;

// Struct to represent the cosine function
struct Cos
{
    static double eval(double x) { return std::cos(x); }
    typedef Neg<Sin> derivative;
};

// Struct to represent the sine function
struct Sin
{
    static double eval(double x) { return std::sin(x); }
    typedef Cos derivative;
};



int main()
{

    typedef ConstInt<4> a; // represents the function a(x) = 4
    typedef ConstInt<5> b; // represents the function b(x) = 5
    typedef Add<a, b> sum; // represents the function sum(x) = a(x) + b(x)
    typedef Mul<a, b> product; // represents the function product(x) = a(x) * b(x)
    typedef Monomial<2> monomial_test; // represents the function monomial(x) = x^2
    typedef Monomial<0> monomial_test2; // represents the function monomial(x) = x^0
    typedef Neg<Monomial<2>> neg_test; // represents the function neg_test(x) = -x^2
    typedef Sin sin_test; // represents the function sin_test(x) = sin(x)
    typedef Cos cos_test; // represents the function cos_test(x) = cos(x)

    // Printing eval() and derivative of the ConstInt object
    std::cout << "Evaluation is: " << ConstInt<5>::eval(0) << ". Its derivative " << ConstInt<5>::derivative::eval(0) << std::endl;

    // Print the value of sum of two ConstInt objects
    std::cout << "Sum of two ConstInt objects is: "<< sum::eval(0) << std::endl;

    // Print the product of two ConstInt objects
    std::cout << "Product of two ConstInt objects is: "<< product::eval(0) << std::endl;

    // Print the evaluation of the monomial function
    std::cout << "Evaluation of the monomial_test function is: " << monomial_test::eval(3) << std::endl;

    // Print the derivative of the monomial function
    std::cout << "Derivative of the monomial_test function is: " << monomial_test::derivative::eval(3) << std::endl;

    // Print the evaluation of the monomial function
    std::cout << "Evaluation of the monomial_test2 function is: " << monomial_test2::eval(3) << std::endl;

    // Print the derivative of the monomial function
    std::cout << "Derivative of the monomial_test2 function is: " << monomial_test2::derivative::eval(3) << std::endl;

    // Testing the Neg struct
    std::cout << "Evaluation of the neg_test function is: " << neg_test::eval(3) << std::endl;
    std::cout << "Derivative of the neg_test function is: " << neg_test::derivative::eval(3) << std::endl;

    // Testing the Sin struct
    std::cout << "Evaluation of the sin_test function is: " << sin_test::eval(0) << std::endl;
    std::cout << "Derivative of the sin_test function is: " << sin_test::derivative::eval(0) << std::endl;

    // Testing the Cos struct
    std::cout << "Evaluation of the cos_test function is: " << cos_test::eval(0) << std::endl;
    std::cout << "Derivative of the cos_test function is: " << cos_test::derivative::eval(0) << std::endl;
    

    return 0;
}

