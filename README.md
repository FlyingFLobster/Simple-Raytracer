![raytraced(21)](https://github.com/FlyingFLobster/Simple-Raytracer/assets/34440763/9dfb0e3d-0d67-4cc1-bb19-11c152d21ddb)

# Simple-Raytracer
 This is a self-learning project in which I am teaching myself the basics of raytracing and its associated algorithms. In its current state this program is only implementing the basics of raytracing with shadows and reflections/refractions and with a single type of geometry, this will be expanded in the future.  
 Currently the locations and properties of all the geometry traced in the output image is determined at the start of the cpp file.  
 There are a number of future features yet to be implemented, such as: Texture support, arbitrary geometry support, more sophisticated bias calculation, brightness attenuation, moveable camera.  
 ## Compilation:
 Requires C++17.  
 g++ raytracing_primer.cpp -o raytracer  
 ## Running:
 In its currently function state, the program is run in the command line and takes in two arguments for the size of the image (width and height as integers), and produces said image in the ppm format.  
