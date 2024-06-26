#include <cmath>
#include <optional>
#include <functional>

// Constants.
#define MAX_RAY_DEPTH 10
#define DEFAULT_INDEX_OF_REFRACTION 1
#define GLASS_INDEX_OF_REFRACTION 1.5

// Structures and Objects.
// Basic 3D Vector struct with inline functions for assignment and vector math operations ala the XNA C# Vector3.
// This is actually based on one of the LEAST efficient forms of a vector in C++, further reading on funky ways to optimise it here: https://www.flipcode.com/archives/Faster_Vector_Math_Using_Templates.shtml
// Operator overloads are so cool.
struct Vector3
{
    float X, Y, Z;

    inline Vector3(void){}

    inline Vector3(const float x, const float y, const float z)
    {X = x; Y = y; Z = z;}

    // Vector3 + Vector3
    inline Vector3 operator + (const Vector3& A) const 
    {return Vector3(X + A.X, Y + A.Y, Z + A.Z);}

    // Vector3 - Vector3
    inline Vector3 operator - (const Vector3& A) const 
    {return Vector3(X - A.X, Y - A.Y, Z - A.Z);}

    // Vector3 Unary -
    inline Vector3 operator - () const 
    {return Vector3(-X, -Y, -Z);}

    // Vector3 + Float (scalar)
    inline Vector3 operator + (const float A) const 
    {return Vector3(X + A, Y + A, Z + A);}

    // Vector3 * Float (scalar)
    inline Vector3 operator * (const float A) const 
    {return Vector3(X * A, Y * A, Z * A);}

    // Vector3 * Int (scalar)
    inline Vector3 operator * (const int A) const 
    {return Vector3(X * A, Y * A, Z * A);}

    // Vector3 . Vector3 (Dot product)
    inline float Dot(const Vector3& A) const
    {return A.X * X + A.Y * Y + A.Z * Z;}

    // Vector3 Magnitude.
    inline float Magnitude() const
    {return sqrt(X*X + Y*Y + Z*Z);}

    // Vector3 Normalisation to unit vector.
    inline Vector3 Normalise() const
    {return Vector3(X / Magnitude(), Y / Magnitude(), Z / Magnitude());} // The compiler shouuuld simplify this upon seeing the repeated expressions.
};

// RGB colour data.
struct Colour
{
    int red, green, blue;

    inline Colour(void){}

    inline Colour(const int in_red, int in_green, int in_blue)
    {red = in_red; green = in_green; blue = in_blue;}

    inline Colour operator + (const Colour& A) const
    {return Colour((red + A.red)/2, (green + A.green)/2, (blue + A.blue)/2);} // Averaging out the values to blend the colours.

    inline Colour operator * (const float& A) const
    {return Colour(static_cast<int>(red * A), static_cast<int>(green * A), static_cast<int>(blue * A));} // Multiply each colour by some scalar.

    // More advanced colour blending function as described on this webpage: https://stackoverflow.com/a/29321264
    Colour BlendColour(Colour A, float ratio)
    {
        int blendRed = sqrt((1 - ratio) * ((red * red) + ratio * (A.red * A.red)));
        int blendGreen = sqrt((1 - ratio) * ((green * green) + ratio * (A.green * A.green)));
        int blendBlue = sqrt((1 - ratio) * ((blue * blue) + ratio * (A.blue * A.blue)));
        return Colour(blendRed, blendGreen, blendBlue);
    }
};

// Object that represents a basic non-specific Ray.
class Ray
{
    private:
        Vector3 m_origin; // The Ray's origin could also be understood as its endpoint.
        Vector3 m_direction; // The Ray's direction is a unit vector and should be normalised on construction.
    public:
        Ray(Vector3 in_origin, Vector3 in_direction)
        {
            m_origin = in_origin;
            m_direction = in_direction.Normalise();
        }

        const Vector3 getOrigin()
        {
            return m_origin;
        }

        const Vector3 getDirection()
        {
            return m_direction;
        }
};

// Function that constructs a Ray specifically for Primary Rays (rays sent from the camera into the scene).
static Ray ConstructPrimaryRay(int rasterX, int rasterY, int width, int height, float aspectRatio, float fov, Vector3 cameraPosition)
{
    // Compute Primary Ray origin and direction.
        //Map the pixel coordinates from Raster to World Space.
        // Pixel (0,0 -> width,height) -> Normalized Device Coordinates(NDC) (0,0 -> 1,1), Keep in mind that in rasterization that NDC is actually [-1,1] in its ranges.
        //float pixelNDCX = (ii + 0.5) / imageWidth; // 0.5 added to both of these so that the we map to the centers of the pixels.
        //float pixelNDCY = (i + 0.5) / imageHeight;
        // NDC -> Screen Space [-1,1], Y coordinates modified so that it is negative for pixels below the X-axis and positive for above.
        //float pixelScreenX = (2 * ndcX) - 1;
        //float pixelScreenY = 1 - (2 * ndcY);
        // Screen Space -> Camera, accounting for aspect ratio and FOV angle.
        //float pixelCameraX = (2 * screenX - 1) * aspectRatio * tan(alpha/2);
        //float pixelCameraY = (1 - (2 * screenY)) * tan(alpha/2);
    // ^Above instructions can be condensed down into just a couple lines like this.
    
    // Careful with the brackets here, messing up the order of operations cost me a few hours of debugging. . .
    float pixelPositionX = (2 * ((rasterX + 0.5) / static_cast<float>(width)) - 1) * aspectRatio * fov;
    float pixelPositionY = (1 - 2 * (rasterY + 0.5) / static_cast<float>(height)) * fov;

    // Final pixel coordinates, keeping in mind that the camera sits exactly one unit away from the image plane on the negative Z-axis.
    Vector3 pixelCameraSpace = {pixelPositionX, pixelPositionY, (cameraPosition.Z - 1)};
    // Vector from the origin of the camera to the pixel on the image plane.
    Vector3 originToPixel = pixelCameraSpace - cameraPosition;
    Vector3 rayDirection = originToPixel.Normalise();

    // Apply camera-to-world transformation matrix here if we want to be able to move the camera around.

    // Finally, we have the Ray's origin and direction.
    return Ray(cameraPosition, rayDirection);
}

// Object that represents a 3D Sphere in the scene.
class Sphere
{
    private:
        bool m_isGlass;
        Colour m_colour;
        Vector3 m_centerPoint;
        float m_radius;
    public:
        Sphere(bool in_isGlass, Colour in_colour, Vector3 in_centerPoint, float in_radius)
        {
            m_isGlass = in_isGlass;
            m_colour = in_colour;
            m_centerPoint = in_centerPoint;
            m_radius = in_radius;
        }

        const bool isGlass()
        {
            return m_isGlass;
        }

        const Colour getColour()
        {
            return m_colour;
        }

        const Vector3 getPosition()
        {
            return m_centerPoint;
        }

        // Returns true if the input ray intersects this sphere, additionally modifies input point references to contain the intersect point and its normal.
        // std::optional<std::reference_wrapper<TYPE>> is used here because normally optional references are not allowed by the standard, but this form is,
        // and by using optionals this way I don't have to further bloat this file with an overload of this decently sized function.
        bool raySphereIntersect(Ray ray, std::optional<std::reference_wrapper<Vector3>> pointHit = std::nullopt, std::optional<std::reference_wrapper<Vector3>> pointNormal = std::nullopt, std::optional<std::reference_wrapper<float>> distance = std::nullopt) 
        {
            // From: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
            // An analytic solution for checking if a Ray intersects a Sphere
            float t_0, t_1; // Possible solutions for if the ray intersects.
            Vector3 L = ray.getOrigin() - m_centerPoint; // Vector from the ray's origin to the center of the sphere.
            
            // Construct the parts of the quadratic equation and then solve it.
            float a = ray.getDirection().Dot(ray.getDirection()); // Ray Magnitude^2.
            float b = 2 * ray.getDirection().Dot(L);
            float c = L.Dot(L) - m_radius * m_radius;

            // Traditional quadratic equation solution ((-b +- sqrt(b^2 - 4ac)) / 2a) can encounter loss of significance errors.
            // This page: https://stackoverflow.com/questions/898076/solve-quadratic-equation-in-c and the above link describe a more robust solution.
            float discriminant = b*b - 4*a*c;
            if (discriminant < 0)
            {
                // If the discriminant is less than 0 then we know there's no (non-complex) solution and can just say there's no intersect.
                return false;
            }
            else if (discriminant == 0)
            {
                // Discriminant == 0 = one solution.
                t_0 = t_1 = -0.5 * b / a; // Simplified since the discriminant being = 0 means there is no reason to add its side of the equation.
            }
            else
            {
                // Discriminant > 0 = two solutions.
                float temp;
                if (b < 0) // Flip the sign depending on b's sign in order to avoid cancellation.
                {
                    temp = -0.5 * (b - sqrt(discriminant));
                }
                else
                {
                    temp = -0.5 * (b + sqrt(discriminant));
                }

                t_0 = temp / a;
                t_1 = c / temp;
            }

            // Swap the intersect points' values if there are two intersections and the first one calculated is further along the ray than the second.
            if (t_0 > t_1)
            {
                std::swap(t_0, t_1);
            }

            // if The first intersection is negative (ie: before the origin of the ray), then discard it for the second intersection.
            if (t_0 < 0)
            {
                if (t_1 < 0) // If both intersections are behind the origin then return false as we shouldn't be able to see it, sort of an error state?
                {
                    return false;
                }
                else
                {
                    t_0 = t_1;
                }
            }

            // Finally, use the first/best intersection solution to calculate the exact point of the intersection and its normal.
            if (pointHit.has_value())
            {
                pointHit.value().get() = (ray.getOrigin() + (ray.getDirection() * t_0));
            }
            if (pointNormal.has_value())
            {
                pointNormal.value().get() = (pointHit.value().get() - m_centerPoint).Normalise();
            }
            if (distance.has_value())
            {
                distance.value().get() = t_0;
            }
            return true;
        }
};

// The class representing light sources in the scene, it is a subclass of sphere to allow for simple construction and added viewability.
class LightSource : public Sphere
{
    private:
        float m_brightness; // Should normally range [0.0, 1.0].
    public:
        LightSource(Colour in_colour, Vector3 in_centerPoint, float in_radius, float in_brightness) : Sphere(false, in_colour, in_centerPoint, in_radius)
        {
            m_brightness = in_brightness;
        }

        float getBrightness()
        {
            return m_brightness;
        }
};

// Function definitions.
Colour Trace(Ray& ray, int depth);