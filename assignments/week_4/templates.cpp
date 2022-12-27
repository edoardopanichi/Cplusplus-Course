#include <iostream>
#include <typeinfo>
#define PRINT_EXPRESSION(expr) std::cout << #expr << ": " << (expr) \
    << " (type: " << typeid(expr).name() << ")" << std::endl

template <typename T> 
T add_simple(const T a, const T b)
{
    return a + b;
}

template <typename T, typename U>
auto add(const T a, const U b)
{
    return a + b;
}

template <typename T>
bool is_int(const T a)
{
    return false;
}

template <>
bool is_int(const int a)
{
    return true;
}

// Number class
template <typename T>
class Number
{
public:
    // Constructor for all types
    Number(const T &value_input)
    : value(value_input) {}

    // Operator overloading
    // SUM
    template <typename U>
    auto operator+(const Number<U> &other) const
    {
        // The type of the return value is the type of the sum of the two values
        return Number<decltype(value + other.value)>(value + other.value);
    }

    // SUBTRACTION
    template <typename U>
    auto operator-(const Number<U> &other) const
    {
        // The type of the return value is the type of the subtraction of the two values
        return Number<decltype(value - other.value)>(value - other.value);
    }

    // MULTIPLICATION
    template <typename U>
    auto operator*(const Number<U> &other) const
    {
        // The type of the return value is the type of the multiplication of the two values
        return Number<decltype(value * other.value)>(value * other.value);
    }

    // DIVISION
    template <typename U>
    auto operator/(const Number<U> &other) const
    {
        // The type of the return value is the type of the division of the two values
        return Number<decltype(value / other.value)>(value / other.value);
    }

    // Attributes
    const T value;
};

// Fibonacci struct
template <int n>
struct fibonacci 
{
  static const int value = fibonacci<n-1>::value + fibonacci<n-2>::value;
};

template <>
struct fibonacci<0> 
{
  static const int value = 0;
};

template <>
struct fibonacci<1> 
{
  static const int value = 1;
};


int main()
{
    int a = 1;
    int b = 2;
    float c = 1.5;
    double d = 2.2;
    double dd = 2.5;
    Number<int> e(1);
    Number<double> f(2.2);

    PRINT_EXPRESSION(add_simple(a, b));
    PRINT_EXPRESSION(add_simple(d, dd));
    PRINT_EXPRESSION(add(a, d));
    PRINT_EXPRESSION(is_int(a));
    PRINT_EXPRESSION(is_int(c));
    
    Number ef = e+f;
    std::cout << "The sum of e+f is: " << ef.value << std::endl;
    Number ef2 = e-f;
    std::cout << "The subtraction of e-f is: " << ef2.value << std::endl;
    Number ef3 = e*f;
    std::cout << "The multiplication of e*f is: " << ef3.value << std::endl;
    Number ef4 = e/f;
    std::cout << "The division of e/f is: " << ef4.value << std::endl;

    std::cout << "The 10th fibonacci number is: " << fibonacci<10>::value << std::endl;
    return 0;
}