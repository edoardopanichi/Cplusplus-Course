import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab

R = np.array([[4, 2, 0, 1],
              [3, 0, 0, 2],
              [0, 1, 1, 1],
              [0, 2, 1, 0]])
A = csc_matrix(R)
b = np.array([-1, -0.5, -1, 2])
x0 = np.array([0, 0, 0, 0])
x, exit_code = bicgstab(A, b, tol=1e-6, maxiter=100, x0=x0)
print(exit_code)
print(x)



# template<typename T>
# // Function to solve the linear system Ax = b using the BiCGStab method
# int bicgstab(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
# {
    
#     // Initialize variables
#     // Initial guess from input: x
#     Vector<T> r_0 = b - A * x; // Initial guessed residual, so the error of the initial guess 
    
#     Vector<T> q_0 = b - A * x; // q_0 is somehow the direction of the search (similar to a gradient)
#     Vector<T> v(maxiter);
#     Vector<T> p(maxiter);
#     v[0] = 0;
#     p[0] = 0;


#     T rho_0 = T(1);
#     T alpha = T(1);
#     T omega_0 = T(1);

#     T rho_k = 0;
#     T beta = 0;
#     T omega_k = 0;
#     Vector<T> h(x.len());
#     Vector<T> s(r_0.len());
#     Vector<T> t(r_0.len());


#     // Iterate
#     for (int i = 0; i < maxiter; i++) 
#     {
#         rho_k = dot(q_0, r_0);
#         beta = (rho_k / rho_0) * (alpha / omega_0);
#         p[i+1] = r_0 + beta * (p[i] - omega_0 * v[i]);
#         v[i+1] = A * p[i+1];
#         alpha = rho_k / dot(q_0, v[i+1]);
#         h = x + alpha * p[i+1];

#         if (norm(b - A * h) < tol) 
#         {
#             x = h;
#             return i;
#         }

#         s = r_0 - alpha * v[i+1];
#         t = A * s;
#         omega_k = dot(t, s) / dot(t, t);
#         x = h + omega_k * s;
#         std::cout << "solution found is: " << std::endl;
#         x.info();

#         if (norm(b - A * x) < tol) {
#             return i;
#         }

#         r_0 = s - omega_k * t;
#         rho_0 = rho_k;
#         omega_0 = omega_k;
#     }

#     // Return number of iterations
#     return -1;
# }



# template<typename T>
# // Function to solve the linear system Ax = b using the BiCGStab method
# int bicgstab(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100)
# {
    
#     // Initialize variables
#     // Initial guess from input: x
#     Vector<T> r_0 = b - A * x; // Initial guessed residual, so the error of the initial guess 
    
#     Vector<T> q_0 = b - A * x; // q_0 is somehow the direction of the search (similar to a gradient)
#     Vector<T> v_0(r_0.len());
#     Vector<T> p_0(r_0.len());
#     T rho_0 = T(1);
#     T alpha = T(1);
#     T omega_0 = T(1);

#     T rho_k = 0;
#     T beta = 0;
#     T omega_k = 0;
#     Vector<T> h(x.len());
#     Vector<T> s(r_0.len());
#     Vector<T> t(r_0.len());


#     // Iterate
#     for (int i = 0; i < maxiter; i++) 
#     {
#         rho_k = dot(q_0, r_0);
#         beta = (rho_k / rho_0) * (alpha / omega_0);
#         p_0 = r_0 + beta * (p_0 - omega_0 * v_0);
#         v_0 = A * p_0;
#         alpha = rho_k / dot(q_0, v_0);
#         h = x + alpha * p_0;

#         if (norm(b - A * h) < tol) 
#         {
#             x = h;
#             return i;
#         }

#         s = r_0 - alpha * v_0;
#         t = A * s;
#         omega_k = dot(t, s) / dot(t, t);
#         x = h + omega_k * s;
#         std::cout << "solution found is: " << std::endl;
#         x.info();

#         if (norm(b - A * x) < tol) {
#             return i;
#         }

#         r_0 = s - omega_k * t;
#         rho_0 = rho_k;
#         omega_0 = omega_k;
#     }

#     // Return number of iterations
#     return -1;
# }
