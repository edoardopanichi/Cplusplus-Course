// Include header file for standard input/output stream library
#include <iostream>

// Your implementation of the swap functions starts here
void swap_by_value(int i1, int i2)
{
    int temp;
    temp = i1;
    i1 = i2;
    i2 = temp; 
}

void swap_by_ref(int& i1, int& i2)
{
    int temp;
    temp = i1;
    i1 = i2;
    i2 = temp; 
}

void swap_by_addr(int* i1, int* i2)
{
    int temp;
    temp = *i1;
    *i1 = *i2;
    *i2 = temp; 
}

// The global main function that is the designated start of the program
int main(int argc, char* argv[]){

    // Initialise integer variables
    int i1=1, i2=2;
    std::cout << "[Initialization] \n" << "value of i1:" << i1 << "  value of i2:" << i2 << "\n" << std::endl;


    // Your testing of the swap functions starts here
    swap_by_value(i1, i2);
    std::cout << "[Swap by value] \n" << "value of i1:" << i1 << "  value of i2:" << i2 << "\n" << std::endl;

    swap_by_ref(i1, i2);
    std::cout << "[Swap by reference] \n" << "value of i1:" << i1 << "  value of i2:" << i2 << "\n" << std::endl;

    swap_by_addr(&i1, &i2);
    std::cout << "[Swap by address] \n" << "value of i1:" << i1 << "  value of i2:" << i2 << "\n" << std::endl;

    // Return code 0 to the operating system (= no error)
    return 0;
}