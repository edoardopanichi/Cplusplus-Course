#include <iostream>
#include <initializer_list>
#include <memory>

class Container
{
public: //[DO NOT MODIFY/REMOVE THESE GETTERS AND SETTERS: THEY ARE USED IN THE SPECTEST]
    int GetLength() const {return length;}
    double* GetData() const {return data;}
    void SetLength(const int length) {this->length = length;}
    void SetData(double* data) {this->data = data;}

public:
    // constructors
    // Default Constructor
    Container() 
    : length(0), data(nullptr)
    { 
        print("\ndefault constructor");
    }
    
    Container(int len)
    : length(len), data(new double[len])
    {
        print("\nConstructor Container(int len)");
    }
    // Constructor from list
    Container(const std::initializer_list<double>& list)
    : length((int)list.size()), data(new double[length])
    {
        std::uninitialized_copy(list.begin(), list.end(), data);
        print("\nConstructor Container(const std::initializer_list<double>& list)");
    }

    //Container(const Container &c) = delete;
    // Copy Constructor
    Container(const Container &other)
    : Container(other.length)
    {
        print("\nCopy Constructor");
        for (int i = 0; i < length; i++)
            data[i] = other.data[i]; 
    }

    // Move constructor
    Container(Container&& other)
    : length(other.length), data(other.data)
    {
        print("\nMove Constructor");
        other.length = 0;
        other.data = nullptr;
    }
    // destructor
    ~Container()
    {
        delete[] data;
        data = nullptr;
        print("\nCalled the destructor");
        length = 0;
    }

    // operators
    void print(const std::string& info) const
    {
        // print the address of this instance, the attributes `length` and
        // `data` and the `info` string
        std::cout << "  " << this << " " << length << " " << data << "  "
            << info << std::endl;
    }

    Container operator=(const Container& other)
    {
        if (this != &other)
        {
            // deleting previous content
            delete[] data;
            // new initialization
            data = new double[other.length];
            length = other.length;
            // copying data
            for (int i = 0; i < length; i++)
                data[i] = other.data[i]; 
        }
        print("\nCopy operator...");
        return *this;
    }

    Container operator=(Container&& other)
    {
        if (this != &other)
        {
            // delete previous content
            delete[] data;

            // Copying the new data. No new initialization of data is needed
            data = other.data;
            length = other.length;

            // Deallocation of other that has been 'consumed' into 'this'. Hence the destructor 
            // does not work on it
            other.data = nullptr;
            other.length = 0;
        }
        print("\nMove operator...");
        return *this;
    }

    Container operator+(const Container &other) const
    {
        Container sum(length);
        for (int i = 0; i < length; i++)
        {
            sum.data[i] = data[i] + other.data[i];
        }
        
        print("\nSum operator...");
        return sum;
    }

private:
    int length;
    double* data;
};

int main()
{

    std::cout << "Container a({ 1, 2, 3 });" << std::endl;
    Container a({ 1, 2, 3 });
    // a({ 1, 2, 3 }) calls the constructor that initialize an object from a list. 
    // Using the list to fill the data attribute.
    std::cout << "  a has address " << &a << std::endl;

    std::cout << "Container b = { 4, 5, 6 };" << std::endl;
    Container b = { 4, 5, 6 };
    // b = { 4, 5, 6 } similarly to the previous one calls the constructor that initialize an 
    // object from a list. If that constructor is made explicit, this definition would not work.
    std::cout << "  b has address " << &b << std::endl;

    std::cout << "Container c(a);" << std::endl;
    Container c(a);
    // c(a) calls the copy constructor that delegates the initialization of c.length to the constructor
    // Container(int length). Then in the copy constructor a.data are copied into c.data.
    std::cout << "  c has address " << &c << std::endl;

    std::cout << "Container d = a + b;" << std::endl;
    Container d = a + b;
    // `a + b`: call `operator+`, which creates a new instance using the
    // constructor `Container(length)`. Then a.data[i]+b.data[i] is assigned to the data 
    // of the new instance.
    std::cout << "  d has address " << &d << std::endl;

    std::cout << "Container e;" << std::endl;
    Container e;
    // The default constructor is called.
    std::cout << "  e has address " << &e << std::endl;
    std::cout << "e = a + b;" << std::endl;
    e = a + b;
    // Same process seen for c = a+b.

    std::cout << "Container f(std::move(a + b));" << std::endl;
    Container f(std::move(a + b));
    // f(std::move(a + b)) calls the move constructor. The move constructor delegates the initialization
    // of f.length to the constructor Container(int length). Then the sum constructor is called for 
    // a+b. The move constructor then copies the data of the sum into f.data. Finally, the destructor
    // of the sum is called.
    std::cout << "  f has address " << &f << std::endl;

    std::cout << "Container g = a + b + c;" << std::endl;
    Container g = a + b + c;
    // g = a + b + c calls the sum operator. The sum operator creates a new instance delegating to the
    // constructor Container(int length). Then the sum operator copies the data of the sum (a+b) into g.data.
    // The process is repeated for the sum (a+b)+c. The destructor of the sum (a+b) is called.
    std::cout << "  g has address " << &g << std::endl;

    std::cout << "Container h;" << std::endl;
    Container h;
    // The default constructor is called.
    std::cout << "  h has address " << &h << std::endl;
    std::cout << "h = a + b + c;" << std::endl;
    h = a + b + c;
    // h = a + b + c works similar to g = a + b + c. But once the final sum is created, the move operator is called.
    // The move operator copies the data of the sum into h.data. Finally, the destructors are called.

    std::cout << "Container i = { a + b + c };" << std::endl;
    Container i = { a + b + c };
    // i = { a + b + c } works as g = a + b + c.
    std::cout << "  i has address " << &i << std::endl;

    std::cout << "end" << std::endl;
    

    // Container x;
    // std::cout << "x has address " << &x << "\n" << std::endl;

    // Container y(12);
    // std::cout << "y length is: " << y.GetLength() << "\n" << std::endl;

    // Container z({1, 3, 4, 5, 7});
    // std::cout << "z length is: " << z.GetLength() << " z content is: " << (z.GetData())[4] << "\n" << std::endl;

    // Container w(z);
    // std::cout << "w length is: " << w.GetLength() << " w content is: " << (w.GetData())[3] << "\n" << std::endl;

    // Container j(std::move(w));
    // std::cout << "j length is: " << j.GetLength() << " j content is: " << (j.GetData())[2] << "\n" << std::endl;

    // Container k = z;
    // std::cout << "k length is: " << k.GetLength() << " k content is: " << (k.GetData())[1] << "\n" << std::endl;

    // Container m = std::move(z);
    // std::cout << "m length is: " << m.GetLength() << " m content is: " << (m.GetData())[1] << "\n" << std::endl;

    // Container n = m + k;
    // std::cout << "n length is: " << n.GetLength() << " n content is: " << (n.GetData())[3] << "\n" << std::endl;

    return 0;
}
