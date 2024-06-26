// #Author Zachary Leon Kidd-Smith
// Heavily based on the teachings and pseudocode provided by www.scratchapixel.com in their lessons on 3d rendering, this work was produced as an exercise
// and may display a lot of weird variations in style (eg: my use of structs with inline functions vs classes in the header file) as I refamiliarise myself with c++.

#include "raytracing_primer.h"
#include <iostream>
#include <cmath>
#include <fstream>

// TODO: Texture support, arbitrary geometry support, more sophisticated bias calculation, brightness attenuation, moveable camera.

// Camera position.
Vector3 cameraPosition = {0, 0, 0};

// Bias sometimes added to the origin of rays to prevent clipping and blinding and other weirdness referred to as "Shadow Acne".
// There are more sophisticated ways with dealing with these phenomena that I can research and add later.
const float bias = 0.0001;

// FOV in degrees, we actually need it in radians but we'll convert it in the preparation phase so that it's easy to modify, could make it adjustable later.
const float fovConst = 90;

// The scene's light source.
LightSource light = LightSource({255, 255, 255}, {0, 20, -20}, 5, 3);

// The scene's background colour.
Colour backgroundColour = Colour{255, 255, 255}; //{50, 50, 80};

// Array containing all spheres in the scene. Sphere = glassY/N, lightsourceY/N, colourRGB, positionV3, Radius.
Sphere spheres[6] = {Sphere(false, Colour{30, 30, 30}, Vector3{0, -315, -40}, 300), // This one acts like a floor because it is so large.
                    Sphere(false, Colour{255, 0, 0}, Vector3{-10, 0, -20}, 3.5),
                    Sphere(false, Colour{0, 0, 230}, Vector3{10, 0, -20}, 3.5),
                    Sphere(false, Colour{0, 120, 160}, Vector3{15, 7, -15}, 1.5),
                    Sphere(true, Colour{0, 90, 0}, Vector3{0, 4, -25}, 3.5),
                    Sphere(false, Colour(0, 255, 0), Vector3(-30, 20, -30), 9)
                    };

int main(int argc, char* argv[])
{
    if (argc <= 2)
    {
        std::cout << "Ray Tracer requires 2 arguments: <Width> and <Height> of output image." << std::endl;
    }
    else
    {
        int imageWidth = atoi(argv[1]);
        int imageHeight = atoi(argv[2]);
        float aspectRatio = static_cast<float>(imageWidth) / imageHeight; // Need to cast int input in order to perform float division instead of int division.
        std::cout << "Raycasting image of dimensions: " << imageWidth << "x" << imageHeight << std::endl;
        std::vector<std::vector<Colour>> pixels(imageHeight, std::vector<Colour>(imageWidth)); // Since the size here is variable based on the user, C++ standard is to use std::vector for arrays, probably safer than C's use of pointers here.
        
        // Calculate the Field Of View (FOV)'s angle (alpha/2) using the triangle formed from the distance between the camera and the viewing plane.
        float fov = tan(fovConst / 2 * M_PI / 180); // PI/180 to convert to radians.

        std::cout << "Begin Ray Trace" << std::endl;
        // Make the height/width loop start here and convert the trace function into a recursive one so that it can easily support multiple reflections/refractions.
        for (int y = 0; y < imageHeight; ++y)
        {
            for (int x = 0; x < imageWidth; ++x)
            {
                // Get direction of the Primary Ray (The Primary Ray being the ray sent from the 'eye' through the pixel (x, y) into the scene).
                Ray primaryRay = ConstructPrimaryRay(x, y, imageWidth, imageHeight, aspectRatio, fov, cameraPosition);
                pixels[y][x] = Trace(primaryRay, 0); // [ii][i] because computer displays are addressed left-to-right top-to-bottom.
            }
        }

        // Write pixels to file directly using a stream in the form of the .ppm format.
        std::ofstream ofs = std::ofstream("./raytraced.ppm", std::ios::out | std::ios::binary);
        ofs << "P6\n" << imageWidth << " " << imageHeight << "\n255\n";
        for (int y = 0; y < imageHeight; ++y)
        {
            for (int x = 0; x < imageWidth; ++x)
            {
                ofs << (unsigned char)pixels[y][x].red << (unsigned char)pixels[y][x].green << (unsigned char)pixels[y][x].blue;
            }
        }
        ofs.close();

    }
}

// Recursive implementation of Ray Tracing Algorithm.
Colour Trace(Ray& ray, int depth)
{
    // Search for Ray-Geometry intersections by looping through each object in the scene and checking if the Ray has intersected them.
    Vector3 pointHit;
    Vector3 normalHit;
    float distance; // The distance between the origin of the ray and the point where it intersects the sphere.
    float minDistance = INFINITY;
    std::optional<Sphere> sphere;
    
    // Real quick check the current depth, doing this first will save a lot of processing if we've already hit the depth limit.
    if (depth > MAX_RAY_DEPTH)
    {
        return backgroundColour;
    }

    for (int i = 0; i < std::size(spheres); ++i)
    {
        if (spheres[i].raySphereIntersect(ray, std::ref(pointHit), std::ref(normalHit), std::ref(distance)))
        {
            if (distance < minDistance)
            {
                // On a successful intersection, update the object and distance so that the above if statement ensures only the closest object is considered hit.
                sphere = spheres[i];
                minDistance = distance;
            }
        }
    }

    // Determine the sphere's material type (diffuse or glass or NULL), return background colour if NULL, 
    // or perform recursive traces into a fresnel equation for reflection/refraction if glass.
    if (!sphere.has_value())
    {
        //std::cout << ++dumbcount << std::endl;
        return backgroundColour;
    }

    if (sphere.value().isGlass() && depth < MAX_RAY_DEPTH)
    {
        // Compute the reflection and refraction rays by recursing back into this function.
        // Would be nice to build constructor functions for these the same way I built one for primary rays, just need to figure out a clean way to separate them.
        Vector3 reflectionDirection = ray.getDirection() - normalHit * 2 * (normalHit.Dot(ray.getDirection())); // From: R = I - 2(N.I)N Where I = Incident vector, N = Normal vector
        Ray reflectionRay = Ray(pointHit + normalHit * bias, reflectionDirection);
        Colour reflectionColour = Trace(reflectionRay, depth + 1);

        // Calculating the refraction direction is a bit more complex, since it changes depending on whether we are entering or exiting the medium (in this case always glass or default(vacuum)).
        // If the ray and the normal of the surface hit are not opposite to eachother (calculated using the dot product), 
        // reverse the normal because we will be inside the sphere for refraction.
        float c_1 = ray.getDirection().Dot(normalHit); // N.I == cos(theta1) == c_1, the angle of incidence.
        float eta_1 = DEFAULT_INDEX_OF_REFRACTION;
        float eta_2 = GLASS_INDEX_OF_REFRACTION;
        if (c_1 < 0)
        {
            // The incident ray is outside the surface, and therefore cos(theta) should to be positive.
            c_1 = -c_1;
        }
        else
        {
            // If the incident ray is coming from inside the surface, then cos(theta) should already be positive but we reverse the normal's direction.
            normalHit = -normalHit;
            // Swap the refraction indices, since we are moving from eta_2 back into eta_1.
            std::swap(eta_1, eta_2);
        }

        float eta = eta_1 / eta_2; // eta_1 / eta_2 == sin(theta2).
        float k = 1 - (eta * eta) * (1 - c_1 * c_1); // Put this in a separate variable from c_2 so that we can use it to check for the critical angle.
        if (k < 0)
        {
            // Angle of the incident is greater than the critical angle, the surface therefore suffers from Total Internal Reflection, 
            //meaning there is no light refracted as it is 100% reflected.
            return reflectionColour;
        }
        else
        { // Otherwise, on with the refraction and fresnel equation.
            float c_2 = sqrt(k); // c_2 == cos(theta2), the angle of refraction.
            Vector3 refractionDirection = ray.getDirection() * eta +  normalHit * (eta * c_1 - c_2); // Final vector form of Snell's Law, pretty cool.
            // Construct and trace the refraction ray.
            Ray refractionRay = Ray(pointHit + normalHit * bias, refractionDirection);
            Colour refractionColour = Trace(refractionRay, depth + 1);

            // Reflection is (obviously) the light reflected, while the refraction is the light transmitted through the new medium.
            // They both make up a ratio depending on the angle of incidence that determines how much of the light either of them use up,
            // always adding up to the amount of light that initially hit the surface since light cannot be destroyed (And this program isn't taking absorption into account so far).
            // Anyway basically the lower the angle of incidence, the greater the amount of light transmitted instead of reflected, and vice versa with reflection and a higher angle of incidence.
            // Further understanding into these two equations for parallel and perpendicular polarised light might require a lot of extra reading, maybe later.
            float fresnelParallel = ((eta_2*c_1 - eta_1*c_2)/(eta_2*c_1 + eta_1*c_2)) * ((eta_2*c_1 - eta_1*c_2)/(eta_2*c_1 + eta_1*c_2));
            float fresnelPerpendicular = ((eta_1*c_2 - eta_2*c_1)/(eta_1*c_2 + eta_2*c_1)) * ((eta_1*c_2 - eta_2*c_1)/(eta_1*c_2 + eta_2*c_1));
            float reflectRatio = (fresnelParallel + fresnelParallel) / 2;
            float refractRatio = 1 - reflectRatio;
            //return reflectionColour * reflectRatio + refractionColour * refractRatio;
            return reflectionColour.BlendColour(refractionColour, refractRatio);
        }
    }
    else if (!sphere.value().isGlass()) // Object is diffuse or opaque, use shadow ray to check for illumination and if it is illuminated return the colour of the object.
    {
        // Send the Shadow Ray from the hit point to the light source to determine the illumination of the point.
        bool isInShadow = false;
        Ray shadowRay = Ray(pointHit + normalHit * bias, (light.getPosition() - pointHit).Normalise());
        for (int i = 0; i < std::size(spheres); ++i)
        {
            if (spheres[i].raySphereIntersect(shadowRay))
            {
                //TODO: ADD CASE FOR IF SPHERE INTERSECTED IS TRANSPARENT.
                isInShadow = true;
                return Colour{0, 0, 0}; // Since we only care if the path to the light source from the point is at all obscured, we can return darkness as soon as we find an object in the way.
            }
        }

        // If the point is being illuminated, return the appropriately lit colour.
        if (!isInShadow)
        {
            // Might need to properly implement Phong lighting here.
            return sphere.value().getColour() * light.getBrightness();
        }
    }
    return backgroundColour; // By default, return the background colour if no intersections are found or max depth is met.
}