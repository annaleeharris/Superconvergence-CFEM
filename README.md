# Superconvergence-CFEM
This folder contains Matlab computer programs to approximate model second order elliptic problems using conforming finte elemet numerical method and a paper explaining about CFEM and superconvergence of CFEM.

The folder contains the folloiwng files:

1. CFEM and Superconvergence CFEM document paper
2. CFEM.m - Matlab computer program computes a numerical approximation u_h.
3. super_CFEM.m - Matlab computer program applies L^2 projection to the existing numerical approximation u_h to enhance the accuracy of u_h and the computational time to calculate u_h.
4. exactu.m - Matlab computer program contains model second order elliptic problems.  It is used to find the difference between the exact solution and the numerical approximation(It calculates the error between u and u_h.)
5. gradientu.m - Matlab computer programs contains model gredients of u.  It is used to find the difference between the exact gradients of u and the numerical approximation gradients u_h (It calculates the error between gradients u and gradients u_h.)
6. rhs.m - Matlab computer program contains f.


