// Include header file for standard input/output stream library
#include <iostream>

// Your implementation of the add functions starts here...
void add(int *array1, int n1, int *array2, int n2, int **array3, int &n3)
{
    // the dimension of the sum array is the max btw n2 and n3
    n3 = std::max(n1, n2);
    int *sum = new int[n3];
    // we initialize array3 such that it points to the pointer that points to the array sum
    *array3 = sum;

    // to accomplish the sum in a unique cycle two variable are defined to switch ON or OFF array1 and array2
    int array1_ON = 1;
    int array2_ON = 1;
    for (int i = 0; i < n3; i++)
    {
        if (i > n1)
        {
            array1_ON = 0;
        }
        if (i > n2)
        {
            array2_ON = 0;
        }
        (*array3)[i] = array1_ON*array1[i] + array2_ON*array2[i];
    }
}

// The global main function that is the designated start of the program
int main(int argc, char *argv[]){

    // Read integer values n1 and n2 from the command line
    int n1, n2;
    if (argc > 2) {
        // If command line arguments are available, use them
        n1 = std::atoi(argv[1]);
        n2 = std::atoi(argv[2]);
    } else {
        // If not, use these default values
        n1 = 10;
        n2 = 15;
    }

    // Allocate and initialize integer arrays array1 and array2
    int *array1 = new int[n1];
    int *array2 = new int[n2];

    int min_n = std::min(n1, n2);
    int max_n = std::max(n1, n2);

    int n3;
    int *array3 = nullptr; 
    
    // up to min_n we fill both arrays.
    for (int i = 0; i < min_n; i++)
    {
        array1[i] = i;
        array2[i] = i;
    }
    // according to which array is longer, we keep filling the longer one.
    if (n1 > n2)
    {
        for (int i = min_n; i < max_n; i++)
        {
            array1[i] = i;
        }
    }
    else
    {
        for (int i = min_n; i < max_n; i++)
        {
            array2[i] = i;
        }
    }
    
    // Test your add function
    add(array1, n1, array2, n2, &array3, n3);
    
    std::cout << "array1: ";
    for (int i = 0; i < n1; i++)
    {
        std::cout << array1[i] << ' ';
    }
    std::cout << "\narray2: ";
    for (int i = 0; i < n2; i++)
    {
        std::cout << array2[i] << ' ';
    }
    std::cout << "\narray3: ";
    for (int i = 0; i < max_n; i++)
    {
        std::cout << array3[i] << ' ';
    }
    std::cout << std::endl;

    //deallocation of the memory;
    delete[] array1;
    delete[] array2;
    delete[] array3;
    array1 = nullptr;
    array2 = nullptr;
    array3 = nullptr;
    // Return code 0 to the operating system (= no error)
    return 0;
}