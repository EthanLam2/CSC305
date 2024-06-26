// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

bool intersect_parallelogram(const Vector3d& ray_origin, const Vector3d& ray_direction,
                              const Vector3d& pgram_origin, const Vector3d& pgram_u, const Vector3d& pgram_v)
{
    MatrixXd D (3,3);
    D << -pgram_u, -pgram_v, ray_direction;
    Vector3d origins = pgram_origin - ray_origin;
    Vector3d lu_decomposition = D.fullPivLu().solve(origins);
    double u = lu_decomposition(0);
    double v = lu_decomposition(1);
    double t = lu_decomposition(2);
    if (u > 1 || u < 0 ) {
        return false;
    } else if (v > 1 || v < 0) {
        return false;
    } else if (t < 0) {
        return false;
    }
    else {
        return true;
    }
    
}

bool intersect_sphere(const Vector3d& ray_origin, const Vector3d& ray_direction,
                      const Vector3d& sphere_center, double sphere_radius)
{
    double a = ray_direction.dot(ray_direction);
    Vector3d oc = ray_origin - sphere_center;
    double b = 2.0 * ray_direction.dot(oc);
    double c = oc.dot(oc) - sphere_radius * sphere_radius;
    double discriminant = b * b - 4 * a * c;    
    if (discriminant < 0) {
        return false;
    }
    double sqrt_discriminant = sqrt(discriminant);
    double t1 = (-b + sqrt_discriminant) / (2 * a);
    double t2 = (-b - sqrt_discriminant) / (2 * a);
    if (t1 < 0 && t2 < 0) {
        return false;
    }   
    return true;
}


void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}


void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            if (intersect_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v))
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection
                MatrixXd D(3,3);
                Vector3d origins = pgram_origin - ray_origin;
                D << -pgram_u, -pgram_v, ray_direction;
                Vector3d lu_decomposition = D.fullPivLu().solve(origins);
                Vector3d ray_intersection = pgram_origin + (lu_decomposition(0) *pgram_u) + (lu_decomposition(1) * pgram_v);
    
                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = (pgram_v.cross(pgram_u)).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);
    

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // TODO: Check if the ray intersects with the parallelogram
            if (intersect_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v))
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection point
                MatrixXd D(3,3);
                Vector3d origins = pgram_origin - ray_origin;
                D << -pgram_u, -pgram_v, ray_direction;
                Vector3d lu_decomposition = D.fullPivLu().solve(origins);
                
                Vector3d ray_intersection = pgram_origin + (lu_decomposition(0) *pgram_u) + (lu_decomposition(1) * pgram_v);

                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = (pgram_v.cross(pgram_u)).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd R = MatrixXd::Zero(800, 800);
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Intersect with the sphere
            // TODO: implement the generic ray sphere intersection
            if (intersect_sphere(ray_origin, ray_direction, sphere_center, sphere_radius))
            {
                // TODO: The ray hit the sphere, compute the exact intersection point
                double a = ray_direction.dot(ray_direction);
                Vector3d oc = ray_origin - sphere_center;
                double b = 2.0 * ray_direction.dot(oc);
                double c = oc.dot(oc) - sphere_radius * sphere_radius;
                double discriminant = (b*b) - (4*a*c);
                double sqrt_discriminant = sqrt(discriminant);
                double t1 = (-b + sqrt_discriminant) / (2 * a);
                double t2 = (-b - sqrt_discriminant) / (2 * a);
                double t = (t1 < t2) ? t1 : t2;
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = ((sphere_center - ray_intersection)/sphere_radius).normalized();

                // TODO: Add shading parameter here
               
                Vector3d l = (ray_intersection - light_position).normalized();
                Vector3d v = (ray_intersection - ray_origin).normalized(); 
                Vector3d h = (v + l).normalized();


                Vector3d diffuse = diffuse_color * std::max(0.0, ray_normal.dot(l));
                Vector3d specular = specular_color * pow(std::max(0.0, ray_normal.dot(h)), specular_exponent);

                // Simple diffuse model
               
                double R_value = ambient + diffuse(0) + specular(0);
                double G_value = ambient + diffuse(1) + specular(1);
                double B_value = ambient + diffuse(2) + specular(2);
                // Clamp to zero
                
                R(i, j) = std::max(R_value, 0.0);
                G(i, j) = std::max(G_value, 0.0);
                B(i, j) = std::max(B_value, 0.0);
                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}