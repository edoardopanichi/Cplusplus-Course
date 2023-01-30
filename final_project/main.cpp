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
        {
            // Setting all the values to zero by default.
            for (int i = 0; i < length; ++i)
            {
                data[i] = 0;
            }
        }
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

    Vector<typename std::common_type<T, U>::type> result(lhs.get_n_rows());
    for (auto it = lhs.get_entries().begin(); it != lhs.get_entries().end(); ++it) 
    {
        int i = it->first.first; // this is the row index of the element of the matrix
        int j = it->first.second; // this is the column index of the element of the matrix
        T value = it->second; // this is the value of the element of the matrix at position (i,j)
        result[i] += value * rhs[j];
    }
    return result;
}


template<typename T>
// Function to solve the linear system Ax = b using the BiCGStab method
int bicgstab(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
{
    
    // Initialize variables

    // VECTORS
    // Initial guess from input: x
    Vector<T> r_0 = b - A * x; // Initial guessed residual, so the error of the initial guess 
    Vector<T> q_0 = r_0; // On wikipedia q_0 is defined as r_hat.

    Vector<T> v(r_0.len());
    Vector<T> p(r_0.len());
    Vector<T> v_0(r_0.len());
    Vector<T> p_0(r_0.len());

    Vector<T> h(x.len()); // h is the first new guess. Each cycle, the guess is updated twice.
    Vector<T> s(r_0.len()); // s is the first new residual. Each cycle, the residual is updated twice.
    Vector<T> t(r_0.len());

    // SCALARS
    T rho_0 = T(1); 
    T alpha = T(1); // alpha is the step size
    T omega_0 = T(1);

    T rho_k = 0;
    T beta = 0;
    T omega_k = 0;

    // Iterate. In each cycle, the guess is updated twice.
    // The first update is saved in h, and the second in x. For the same reason, the residual is updated twice.
    for (int i = 0; i < maxiter; ++i) 
    {
        rho_k = dot(q_0, r_0);
        beta = (rho_k / rho_0) * (alpha / omega_0);
        p = r_0 + beta * (p_0 - omega_0 * v_0); // p is the direction of the search
        v = A * p;
        alpha = rho_k / dot(q_0, v); // update of the step size
        h = x + alpha * p; // update of the guess
        
        if (norm(b - A * h) < tol)
        {
            x = h;
            // Print solution before exiting
            std::cout << "Solution Found is: " << std::endl;
            x.info();
            return i;
        }

        s = r_0 - alpha * v; // s is the new residual
        t = A * s; 
        omega_k = dot(t, s) / dot(t, t);
        x = h + omega_k * s; // update of the guess

        if (norm(b - A * x) < tol) 
        {
            // Print solution before exiting
            std::cout << "Solution Found is: " << std::endl;
            x.info();
            return i;
        }

        // new values become old values for the next iteration
        r_0 = s - omega_k * t; 
        p_0 = p;
        v_0 = v;
        rho_0 = rho_k;
        omega_0 = omega_k;        
    }

    // Return number of iterations
    return -1;
}

// Function to implement the Heun integration method
template<typename T>
// f is the vector of functions to integrate, y is the vector of initial conditions, h is the step size, and t is the initial time
void heun(const Vector<std::function<T(Vector<T> const&, T)> >& f, Vector<T>& y, T h, T& t)
{
    // The method is divided in two steps. The first step is to compute the slope at the current point, like in Euler's method.
    // The second step is to compute the slope at the point where the first step ends, and then compute the average of the two slopes.
    // The average is then used to update the current point.
    
    // FIRST STEP
    Vector<T> k_1 = f(y, t); // slope at the current point
    Vector<T> y_1 = y + h * k_1; // point where the first step ends (Euler's method)

    // SECOND STEP
    Vector<T> k_2 = f(y_1, t + h); // slope at the point where the first step ends
    Vector<T> k_avg = (k_1 + k_2) / 2; // average of the two slopes
    y = y + h * k_avg; // update of the current point
    t += h; // update of the current time
    
};

// Class to implement a model of a 2-DoF Walker. 
template<typename T>
class SimplestWalker
{
    public:
        // CONSTRUCTORS
        // Constructor with initial condition (y0 = y(t0)), initial time and slope of the path
        SimplestWalker(const& Vector<T> y0, T t0, T gamma) 
        :   y(y0), 
            t(t0), 
            slope(gamma)
        {}

        // METHODS
        // Function to compute the derivative of the position
        Vector<T> derivative(const Vector<T>& y) const
        {
            Vector<T> y_dot(y.len()); // it contains the derivative of y.
            
            // Value of the derivative computed on paper starting from the equations of motion
            y_dot[0] = y[3]; 
            y_dot[1] = y[4];
            y_dot[2] = sin(y[2]-slope) + y[4]^2 * sin(y[1]) - cos(y[2]-slope) * sin(y[1]);
            y_dot[3] = sin(y[2]-slope);

            return y_dot;
        }

        // Function to compute the evolution of the system after a time step h
        const Vector<T>& step(T h)
        {
            // The system is integrated using the Heun method
            heun(derivative, y, h, t);
            return y;
        }


    private:
        // ATTRIBUTES
        Vector<T> y; // Vector of the current position
        T t; // Current time
        T slope; // Slope of the path
    public:
        Matrix<T> A; // Matrix of the system
        Vector<T> b; // Vector of the right hand side of the system
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

    // // Testing the get function of matrix class
    // std::cout << "B.get(0, 0) = " << B.get({0, 0}) << std::endl;

    // Testing the * operator of matrix class with a vector
    Vector<double> v14({1, 2, 3});
    auto v15 = B * v14;
    v15.info();

    // Testing the bicgstab function
    Matrix<double> A2(4, 4);
    A2[{0, 0}] = 2.0;
    A2[{0, 1}] = 2.0;
    // A2[{0, 2}] = 0.0;
    A2[{0, 3}] = 1.0;

    A2[{1, 0}] = 3.0;
    // A2[{1, 1}] = 4.0;
    // A2[{1, 2}] = 1.0;
    A2[{1, 3}] = 2.0;

    // A2[{2, 0}] = 0.0;
    A2[{2, 1}] = 1.0;
    A2[{2, 2}] = 1.0;
    A2[{2, 3}] = 1.0;
    
    // A2[{3, 0}] = 0.0;
    A2[{3, 1}] = 2.0;
    A2[{3, 2}] = 1.0;
    // A2[{3, 3}] = 4.0;
    A2.info();
    Vector<double> b2({-1, -0.5, -1, 3});
    Vector<double> x2({0, 0, 0, 0});

    auto iterations = bicgstab(A2, b2, x2, 1e-6);
    std::cout << "bicgstab iterations: " << iterations << std::endl;
    
    // Your testing of the simplest walker class starts here

    return 0;
}