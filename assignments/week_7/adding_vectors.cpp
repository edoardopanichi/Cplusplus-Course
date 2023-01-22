#ifdef _MSC_VER
#  define NOINLINE __declspec(noinline)
#else
#  define NOINLINE __attribute__((noinline))
#endif

#include <chrono>
#include <vector>

// Function to add the i-th element of each vector in the argument pack
template<typename V, typename... T>
V add_elems(const int i, const std::vector<V>& head, const T&... tail)
{
    if constexpr (sizeof...(tail) > 0)
        return head[i] + add_elems(i, tail...);
    else
        return head[i];
}

// TODO: Replace the following function template. You may change the template
// arguments and function arguments if necessary.
// add_vectors_singleloop should be implemented using a single loop and exploiting add_elems.
template<typename V, typename... T>
std::vector<V> add_vectors_singleloop(const std::vector<V>& head, const T&... tail)
{
    std::vector<V> result(head.size());
    for (unsigned long int i = 0; i < head.size(); i++)
        result[i] = add_elems(i, head, tail...);
    return result;
}

// TODO: Replace the following function template. You may change the template
// arguments and function arguments if necessary.
template<typename V, typename... T>
std::vector<V> add_vectors_simple(const std::vector<V>& head, const T&... tail)
{
    std::vector<V> result(head.size());
    if constexpr (sizeof...(tail) > 0)
    {
        // Recursively call add_vectors_simple on the tail of the argument pack
        auto tail_result = add_vectors_simple(tail...);
        // The for loop below adds the head and tail for each element of the vector
        for (unsigned long int i = 0; i < head.size(); i++)
            result[i] = head[i] + tail_result[i];
    }
    else
    {
        for (unsigned long int i = 0; i < head.size(); i++)
            result[i] = head[i];
    }
    return result;
}



NOINLINE std::vector<double> test_add_vectors_singleloop(
    const std::vector<double>& a, const std::vector<double>& b,
    const std::vector<double>& c, const std::vector<double>& d)
{
    return add_vectors_singleloop(a, b, c, d);
}

NOINLINE std::vector<double> test_add_vectors_simple(
    const std::vector<double>& a, const std::vector<double>& b,
    const std::vector<double>& c, const std::vector<double>& d)
{
    return add_vectors_simple(a, b, c, d);
}

#include <iostream>
#include <cstring>

int main(int argc, char* argv[])
{
    int n = 10000;
    if (argc > 1) n = std::atoi(argv[1]);

    std::vector<double> a(n);
    std::vector<double> b(n);
    std::vector<double> c(n);
    std::vector<double> d(n);

    {
        std::cout << "Testing simple version........";
        auto tic = std::chrono::system_clock::now();
        for (int i = 0; i < 100; i++)
            test_add_vectors_simple(a, b, c, d);
        auto toc = std::chrono::system_clock::now();
        std::cout << std::chrono::duration<double>(toc-tic).count() << "s" << std::endl;
    }

    {
        std::cout << "Testing single-loop version...";
        auto tic = std::chrono::system_clock::now();
        for (int i = 0; i < 100; i++)
            test_add_vectors_singleloop(a, b, c, d);
        auto toc = std::chrono::system_clock::now();
        std::cout << std::chrono::duration<double>(toc-tic).count() << "s" << std::endl;
    }

    double sum_elem = add_elems(0, a, b, c, d);
    std::cout << "Sum of first element: " << sum_elem << std::endl;
    return 0;
}