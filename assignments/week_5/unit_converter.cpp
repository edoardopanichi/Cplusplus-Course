// In this file the Measure_add trait struct is implemented in this way:
// The unit with the biggest scale is converted to the unit with the smallest scale and then 
// the two values are added. The result is in the unit with the smallest scale.

#include <iostream>

// Implementation of the enum class Unit
enum class Unit { km, m, cm };

// Measure trait struct
template <int v, Unit u>
struct Measure
{
    static const int value = v;
    static const Unit unit = u;
};


// Implementing a trait the computes the scaling factor between two units
template <Unit T, Unit U>
struct scale_factor
{
    static const int value_1 = 1;
    static const int value_2 = 1;
    static const Unit unit = T;
};

template <>
struct scale_factor<Unit::km, Unit::m>
{
    static const int value_1 = 1000;
    static const int value_2 = 1;
    static const Unit unit = Unit::m;
};

template <>
struct scale_factor<Unit::km, Unit::cm>
{
    static const int value_1 = 100000;
    static const int value_2 = 1;
    static const Unit unit = Unit::cm;
};

template <>
struct scale_factor<Unit::m, Unit::km>
{
    static const int value_1 = 1;
    static const int value_2 = 1000;
    static const Unit unit = Unit::m;
};

template <>
struct scale_factor<Unit::m, Unit::cm>
{
    static const int value_1 = 100;
    static const int value_2 = 1;
    static const Unit unit = Unit::cm;
};

template <>
struct scale_factor<Unit::cm, Unit::km>
{
    static const int value_1 = 1;
    static const int value_2 = 100000;
    static const Unit unit = Unit::cm;
};

template <>
struct scale_factor<Unit::cm, Unit::m>
{
    static const int value_1 = 1;
    static const int value_2 = 100;
    static const Unit unit = Unit::cm;
};

// Measure_add trait struct: converts to the smallest unit and adds
template <typename T, typename U>
struct Measure_add
{   
    static const int value = (T::value * scale_factor<T::unit, U::unit>::value_1) + (U::value * scale_factor<T::unit, U::unit>::value_2);
    static const Unit unit = scale_factor<T::unit, U::unit>::unit;
};


int main()
{
    std::cout << Measure_add< Measure<1,Unit::km>, Measure<1,Unit::m> >::value  << std::endl;
    std::cout << Measure_add< Measure<10,Unit::m>, Measure<20,Unit::cm> >::value << std::endl;

    return 0;
}