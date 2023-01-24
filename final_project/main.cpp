#include <cmath>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

template <typename T>
class Vector
{
    public:
        // Default constructor that sets the length to zero
        Vector() 
        :   length(0),
            data(nullptr)
        {} 
        // Constructor that takes a length as input
        Vector(const int& length)
        :   length(length),
            data(new T[length])
        {}
        // Copy constructor
        Vector(const Vector& other)
        :   length(other.length),
            data(new T[other.length])
        {
            for (int i = 0; i < length; ++i)
            {
                data[i] = other.data[i];
            }
        }
        // Move constructor
        Vector(Vector&& other)
        :   length(other.length),
            data(other.data)
        {
            other.length = 0;
            other.data = nullptr;
        }
        // Initializer list constructor
        Vector(const std::initializer_list<T>& list)
        :   Vector((int)list.size())
        {
            std::uninitialized_copy(list.begin(), list.end(), data);
        } 

        // Operator overloading

        // Copy assignment operator
        Vector& operator=(const Vector& other)
        {
            if (this != &other)
            {
                delete[] data;
                length = other.length;
                data = new T[length];
                for (int i = 0; i < length; ++i)
                {
                    data[i] = other.data[i];
                }
            }
            return *this;
        }
        // Move assignment operator
        Vector& operator=(Vector&& other)
        {
            if (this != &other)
            {
                delete[] data;
                length = other.length;
                data = other.data;
                other.length = 0;
                other.data = nullptr;
            }
            return *this;
        }
        // [] operator
        T& operator[](int i)
        {
            if (i < 0 || i >= length)
            {
                throw std::out_of_range("Index out of range");
            }
            return data[i];
        }
        // [] operator (const)
        const T& operator[](int i) const
        {
            if (i < 0 || i >= length)
            {
                throw std::out_of_range("Index out of range");
            }
            return data[i];
        }
        // + operator that adds two vectors of any type and returns a vector of the common type
        template<typename U>
        auto operator+(const Vector<U>& other) const
        {
            if (length != other.length)
            {
                throw std::invalid_argument("Vectors must have the same length");
            }
            Vector<typename std::common_type<T,U>::type> result(length);
            for (int i = 0; i < length; ++i)
            {
                result[i] = data[i] + other[i];
            }
            return result;
        }
        // - operator that subtracts two vectors of any type and returns a vector of the common type
        template<typename U>
        auto operator-(const Vector<U>& other) const
        {
            if (length != other.length)
            {
                throw std::invalid_argument("Vectors must have the same length");
            }
            Vector<typename std::common_type<T,U>::type> result(length);
            for (int i = 0; i < length; ++i)
            {
                result[i] = data[i] - other[i];
            }
            return result;
        }
        // * operator that multiplies a vector with a scalar of any type and returns a vector of the common type
        template<typename U>
        auto operator*(const U& scalar) const
        {
            Vector<typename std::common_type<T,U>::type> result(length);
            for (int i = 0; i < length; ++i)
            {
                result[i] = data[i] * scalar;
            }
            return result;
        }
        // * operator that multiplies scalar and vector, and checks that scalar is arithmetic
        template<typename U, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>
        friend auto operator*(const U& scalar, const Vector<T>& vec)
        {
            return vec * scalar;
        }

        // Destructor
        ~Vector()
        {
            delete[] data;
            length = 0;
        }

        // Attributes
        int length;
        T* data;

        // Methods
        void info() const
        {
            std::cout << "Length: " << length << "  ";
            std::cout << "Data: ";
            for (int i = 0; i < length; ++i)
            {
                std::cout << data[i] << " ";
            }
            std::cout << std::endl;
        }

        // To get the length of the vector
        const int len() const
        {
            return length;
        }
};

template<typename T, typename U>
typename std::common_type<T,U>::type 
dot(const Vector<T>& lhs, 
    const Vector<U>& rhs)
{
    if (lhs.length != rhs.length)
    {
        throw std::invalid_argument("Vectors must have the same length");
    }
    typename std::common_type<T,U>::type result = 0;
    for (int i = 0; i < lhs.length; ++i)
    {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template<typename T>
T norm(const Vector<T>& vec)
{
    // Function that returns the norm of a vector
    return std::sqrt(dot(vec, vec));
}

// Class to represent a sparse matrices using a map
template <typename T>
class Matrix
{
    public:
        //Constructors
        Matrix()
        :   n_rows(0), n_cols(0)
        {}

        Matrix(int n_rows, int n_cols)
        :   n_rows(n_rows), n_cols(n_cols)
        {}

        // Destructor
        ~Matrix()
        {
            n_rows = 0;
            n_cols = 0;
        }

        // Operator Overloading
        // [] operator
        T& operator[](const std::pair<int, int>& ij)
        {
            if (ij.first < 0 || ij.first >= n_rows || ij.second < 0 || ij.second >= n_cols)
            {
                throw std::out_of_range("Index out of range");
            }
            return entries[ij];
        }

        // () operator, returns the value of the entry ij as a const reference
        const T& operator()(const std::pair<int, int>& ij) const
        {
            if (ij.first < 0 || ij.first >= n_rows || ij.second < 0 || ij.second >= n_cols)
            {
                throw std::out_of_range("Index out of range");
            }
            return entries.at(ij);
        }

        // Methods
        // Print all the information of the matrix
        void info() const
        {
            std::cout << "Number of rows: " << n_rows << std::endl;
            std::cout << "Number of columns: " << n_cols << std::endl;
            std::cout << "Entries: " << std::endl;
            for (auto it = entries.begin(); it != entries.end(); ++it)
            {
                std::cout << "(" << it->first.first << ", " << it->first.second << ") = " << it->second << std::endl;
            }
        }

        // Function to get an element of the matrix, using an iterator over all the entries of the map
        // const T get(const std::pair<int, int>& ij) const
        // {
        //     if (ij.first < 0 || ij.first >= n_rows || ij.second < 0 || ij.second >= n_cols)
        //     {
        //         throw std::out_of_range("Index out of range");
        //     }
        //     for (auto it = entries.begin(); it != entries.end(); ++it) 
        //     {
        //         int i = it->first.first;
        //         int j = it->first.second;
        //         found_value = it->second;
        //         if (i == ij.first && j == ij.second) 
        //             // exit the loop if the element is found
        //             return found_value;    
        //     }
        //     // return a default-constructed T if the element is not found
        //     return T();
        // }

        // Function to get the number of rows
        const int get_n_rows() const
        {
            return n_rows;
        }

        // Function to get the number of columns
        const int get_n_cols() const
        {
            return n_cols;
        }
        
        // Function to get the entries
        const std::map<std::pair<int, int>, T>& get_entries() const
        {
            return entries;
        }

    private:
        // Attributes
        std::map<std::pair<int, int>, T> entries;
        int n_rows;
        int n_cols;
        mutable T found_value;

};

template<typename T, typename U>
Vector<typename std::common_type<T,U>::type> operator*(const Matrix<T>& lhs, const Vector<U>& rhs)
{
    
    if (lhs.get_n_cols() != rhs.length)
    {
        throw std::invalid_argument("Matrix and vector must have compatible dimensions");
    }
    Vector<typename std::common_type<T,U>::type> result(lhs.get_n_rows());
    for (auto it = lhs.get_entries().begin(); it != lhs.get_entries().end(); ++it) 
    {
        //std::cout << *it << std::endl;
        int i = it->first.first;
        int j = it->first.second;
        T value = it->second;
        result[i] += value * rhs[j];
    }
    return result;
}


template<typename T>
int bicgstab(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
{
    // Your implementation of the bicgstab function starts here
    return 0;
}

template<typename T>
void heun(const Vector<std::function<T(Vector<T> const&, T)> >& f, Vector<T>& y, T h, T& t)
{
    // Your implementation of the heun function starts here
};

template<typename T>
class SimplestWalker
{
    // Your implementation of the simplest walker class starts here
};

int main(int argc, char* argv[])
{   
    Vector<double> v1; // Default constructor
    Vector<double> v2(3); // Constructor with length
    Vector<double> v3(4); // Vector to be moved into v4
    Vector<double> v4(std::move(v3)); // Move constructor
    Vector<double> v5(v2); // Copy constructor
    Vector<double> v6({1, 2, 3, 4, 5}); // Initializer list constructor

    // Printing the content of the vectors.
    std::cout << "v1: " << std::endl;
    v1.info();
    std::cout << "v2: " << std::endl;
    v2.info();
    std::cout << "v3: " << std::endl;
    v3.info();
    std::cout << "v4: " << std::endl;
    v4.info();
    std::cout << "v5: " << std::endl;
    v5.info();
    std::cout << "v6: " << std::endl;
    v6.info();

    v5 = v6; // Copy assignment operator

    std::cout << "v5 after copy assignment: " << std::endl;
    v5.info();

    v6 = std::move(v5); // Move assignment operator

    std::cout << "v5 after move assignment: " << std::endl;
    v5.info();
    std::cout << "v6 after move assignment: " << std::endl;
    v6.info();

    // Testing the + operator
    Vector<double> v7({1, 2, 3, 4, 5});
    Vector<int> v8({1, 2, 3, 4, 5});
    auto v9 = v7 + v8 - v7;
    std::cout << "v9: " << std::endl;
    v9.info();

    // Testing the * operator
    auto v10 = v9 * 2.0;
    std::cout << "v10: " << std::endl;
    v10.info();

    auto v11 = 2.0 * v9;
    std::cout << "v11: " << std::endl;
    v11.info();
    std::cout << "v11 len: " << v11.len() << std::endl;

    // Testing the dot function
    auto v12 = Vector<double>({1, 2, 3});
    auto v13 = Vector<double>({3, 1, 2});
    auto dot_result = dot(v12, v13);
    std::cout << "dot_result: " << dot_result << std::endl;

    // Testing the norm function
    auto norm_result = norm(v12);
    std::cout << "norm_result: " << norm_result << std::endl;

    // Testing the class Matrix
    Matrix<double> A;
    Matrix<double> B(3, 3);
    A.info();
    B[{0, 0}] = 1.0;
    B[{0, 1}] = 2.0;
    B[{0, 2}] = 3.0;
    B[{1, 0}] = 4.0;
    B[{1, 1}] = 5.0;
    B[{1, 2}] = 6.0;
    B[{2, 0}] = 7.0;
    B[{2, 1}] = 8.0;
    B[{2, 2}] = 9.0;

    B.info();

    // Testing () operator of matrix class
    std::cout << "B(0, 0) = " << B({0, 0}) << std::endl;

    // Testing the get function of matrix class
    std::cout << "B.get(0, 0) = " << B.get({0, 0}) << std::endl;

    // Testing the * operator of matrix class with a vector
    Vector<double> v14({1, 2, 3});
    auto v15 = B * v14;
    v15.info();

    
    // Your testing of the simplest walker class starts here

    return 0;
}