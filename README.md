# RotorIy-Standalone : Rotor Identification with IMU using geometric algebra for standalone envilonments

The original RotoIy program uses Gaigen2.5 GA implementation generator to generate GA computaion. It's fast but too heavy for standalone environments with microcontrollers.

The current version is highly experimental.

The associated main.cpp is only a sample and isn't intended to the real one, though it's tested with hachidori/bee3b WiFi remote sensor. The identification is done with RotorIy files.

The rotor identification is based on the method proposed at

[Rotation Identification in Geometric Algebra: Theory and Application to the Navigation of Underwater Robots in the Field](https://pdfs.semanticscholar.org/b2a3/a6b7221b215840a7f910179665c8419b0ec0.pdf)

When fuising magnetometer too, Fontijne-Dost algorism of 3D rotor reconstruction is used to get rotor P which maps x-axis to the north and z-axis to the opposite direction of the gravity.

GA computations are compiled to vecotor/bivector arithmetics by hand, augumented with Mathematica. There are 2 implementations in RotorIy. The class RotorVS is based versor(rotor) multiplication and RotorBV is based the bivector updating method which is given Candy and Lasenby.

Some trigonometric functions - sinf, conf, atan2f and sincf functions are re-implemented for the limited standalone environments. sinf is computed with the table lookup and a few multiplications.

