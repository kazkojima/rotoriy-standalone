# RotorIy-Standalone : Rotor Identification with IMU using geometric algebra for standalone environments

The original RotoIy program uses Gaigen2.5 GA implementation generator to generate GA computation. It's fast but too heavy for standalone environments with micro-controllers.

The current version is highly experimental.

The associated main.cpp is only a sample and isn't intended to the real one, though it's tested with hachidori/bee3b WiFi remote sensor. The identification is done with RotorIy files.

The rotor identification is based on the method proposed at

[Rotation Identification in Geometric Algebra: Theory and Application to the Navigation of Underwater Robots in the Field](https://pdfs.semanticscholar.org/b2a3/a6b7221b215840a7f910179665c8419b0ec0.pdf)

When fusing magnetometer too, Fontijne-Dost algorithm[¹] of 3D rotor reconstruction is used to get rotor P which maps x-axis to the north and z-axis to the opposite direction of the gravity.

GA computations are compiled to vector/bivector arithmetic's by hand, augmented with Mathematica. There are 2 implementations in RotorIy. The class RotorVS is based versor(rotor) multiplication and RotorBV is based the bivector updating method which is given Candy and Lasenby.[²]

Some trigonometric functions - sinf, conf, atan2f and sincf functions are re-implemented for the limited standalone environments. sinf is computed with the table lookup and a few multiplications.

[¹] Daniel Fontijne, Leo Dorst: Reconstructing Rotations and Rigid Body Motions from Exact Point Correspondences Through Reflections. Guide to Geometric Algebra in Practice 2011: 63-78

[²] Liam Candy, Joan Lasenby: Attitude and Position Tracking. Guide to Geometric Algebra in Practice 2011: 105-125
