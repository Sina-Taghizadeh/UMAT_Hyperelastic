# UMAT_Hyperelastic
Comparison of Abaqus Library, UMAT, and Analytical Solution for Neo-Hookean Hyperelastic Material Under Large Deformations
# Problem Definition
To enable analytical solution, we need to investigate problems where the deformation gradient tensor is known precisely. We can then substitute it into the formula for stress as a function of the deformation gradient and obtain the stress tensor, which is the same for the entire domain except at the support points in these cases. For this purpose, we choose two important problems in solid mechanics: simple shear and simple extension.
In both loading cases, the problem we consider is a beam with a length of 100 mm and a cross-sectional area of 10*10 mm2, made of a Neo-Hookean hyperelastic material with a shear modulus of 1 MPa and a bulk modulus of 0.78 MPa. We also choose a shear strain of 0.4 for the simple shear case and a longitudinal strain of 0.2 for the simple extension case to ensure that we are well within the large deformation regime. Therefore, the deformation gradient tensors for these cases can be represented as the following matrices:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/44c52ada-9f91-4e18-b349-bd216daa06cd)

We obtain the response for these two problems using three methods:
# 1. Exact Analytical Solution
We have various constitutive relations for Neo-Hookean materials. The constitutive relation used by Abaqus software, according to its help section, is as follows:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/ec20dafe-5c9e-47cc-815c-c334d52fa304)

where J is the determinant of the deformation gradient tensor, and Bbar, which is called the deviatoric left Cauchy-Green deformation tensor, is defined as follows:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/aaac1a16-988d-46b3-9c8f-c9bf3b264eea)

where Fbar, the distortion gradient, is defined as follows:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/fa95f2ba-3c26-45ff-8fc4-99ef187c76e3)

C10 and D1 are also calculated based on the shear modulus μ0 and the bulk modulus K0 using the following relations:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/f892d19f-9d25-46e6-b8af-2c33fe34fcf4)

Thus, for the two aforementioned problems, we can obtain their stress tensors by substituting their deformation gradients into the above relations. For simplicity, this formula has been implemented in MATLAB and is located in the file 'AnalyticalNeoHookean.m' in the repository. The resulting stress tensors are as follows:

![image](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic/assets/162900845/b0495973-75b3-432d-9a4c-023533a08b04)

# 2. Abaqus Library Solution
In this section, we solve the problem by modeling it and using the Neo-Hookean hyperelastic material available in the Abaqus Complete Abaqus Environment (CAE), which is pre-written by its designers. Since we are dealing with large deformations, we turn on the "nonlinear geometry" option in the step section and reduce the "initial increment" value to a much smaller value to start with a small value and automatically increase it if needed. We implement the two mentioned loading conditions in the "Load" section. The inp files for Simple Extension and Simple Shear cases are named 'NeoHookeanSimpleExtensionCAE.inp' and 'NeoHookeanSimpleShearCAE.inp', respectively, and are located in the repository. After performing the analysis, we see that the resulting stress tensor is "completely consistent" with the analytical solution.

# 3. Solution using UMAT subroutine
Although it is recommended to use the UHYPER subroutine for hyperelastic materials and define the strain energy density and its derivatives in it, here I used the UMAT subroutine and defined the Cauchy stress tensor and the system Jacobian matrix, which is known as DDSDDE in this environment. The subroutine written in Fortran is placed in the file named 'CompresibleNeoHookean.for' in the repository. The inp files for the Simple Extension and Simple Shear cases, which use the written UMAT subroutine, are named 'NeoHookeanSimpleExtensionUMAT.inp' and 'NeoHookeanSimpleShearUMAT.inp', respectively, and are placed in the repository. After performing the analysis in this case, we also see that the resulting stress tensor is "completely consistent" with the analytical solution and the CAE solution.

# Acknowledgments
A large part of the Fortran code written for the Neo-Hookean material was inspired by the code written by Dr. Milad Vahedian, which can be found at https://mecademy.org/abaqus-subroutine/. For more information on this subroutine and other subroutines, you can refer to his free tutorials (in Persian) on this site.
