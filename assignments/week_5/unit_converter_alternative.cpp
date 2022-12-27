// In this file the Measure_add trait struct is implemented in this way:
// Both units are converted into cm and then added. The result is in cm.

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

// Measure_add trait struct: converts to the centimeters and adds
template <typename T, typename U>
struct Measure_add
{
    // In this case, the T::unit == Unit::km condition is checked first. If it is true, 
    // then 100000 is returned as the value of T_converted. If it is false, the 
    // T::unit == Unit::m condition is checked next. If this condition is true, 100 is returned 
    // as the value of T_converted. If both of these conditions are false, then 1 is returned as 
    // the value of T_converted.
    static const int T_converted = T::unit == Unit::km ? 100000 : T::unit == Unit::m ? 100 : 1;
    // T_converted and U_converted have to be static const int, because they are used in static
    // context (in the definition of the value member of the Measure_add struct). Basically,
    // when if T_converted and U_converted are not static const int, the compiler will complain
    // that the value of T_converted and U_converted is not known at compile time. Indeed these two
    // variables, not static, exists only when an object of the Measure_add struct is created.
    // Whereas static const int variables are known at compile time.
    static const int U_converted = U::unit == Unit::km ? 100000 : U::unit == Unit::m ? 100 : 1;
    static const int value = (T::value * T_converted) + (U::value * U_converted);
};


int main()
{
    std::cout << Measure_add< Measure<10,Unit::m>, Measure<20,Unit::m> >::value  << std::endl;
    std::cout << Measure_add< Measure<10,Unit::m>, Measure<20,Unit::cm> >::value << std::endl;

    return 0;
}