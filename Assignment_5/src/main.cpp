// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5;       //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    // TODO: Setup uniform

    // TODO: Setup camera, compute w, u, v
    const Vector3d w = -camera_gaze.normalized();
    const Vector3d u = camera_top.cross(w).normalized();
    const Vector3d v = w.cross(u);
    // TODO: Compute the camera transformation
    Matrix4d matrixCamera;

    matrixCamera << u(0), v(0), w(0), camera_position(0),
                    u(1), v(1), w(1), camera_position(1),
                    u(2), v(2), w(2), camera_position(2),
                    0, 0, 0, 1;


    Matrix4d transformedCamera;
    transformedCamera = matrixCamera.inverse();

    // TODO: Setup projection matrix

    const double t = near_plane * tan(field_of_view / 2.0);
    const double r = t * aspect_ratio;
    const double l = -r;
    const double b = -t;
    const double n = -near_plane;
    const double f = -far_plane;

    Matrix4d orthMatrix;
    orthMatrix << 2 / (r - l), 0, 0, -(r + l) / (r - l), 
                  0, 2 / (t - b), 0, -(t + b) / (t - b), 
                  0, 0, 2 / (n - f), -(n + f) / (n - f), 
                  0, 0, 0, 1;

    Matrix4d perspectiveCam;
    if (is_perspective)
    {
        // TODO setup perspective camera
        perspectiveCam << n, 0, 0, 0,
                          0, n, 0, 0,
                          0, 0, n + f, -f * n,
                          0, 0, 1, 0;

        uniform.cam = orthMatrix * perspectiveCam * transformedCamera;
    }
    else
    {
        uniform.cam = orthMatrix * transformedCamera;
    }
}
Matrix4d compute_rotation(const double alpha)
{
    //TODO: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    Matrix4d res;
    double sinAlpha = sin(alpha);
    double cosAlpha = cos(alpha);
    
    res << cosAlpha, 0, sinAlpha, 0,
    0, 1, 0, 0,
    -sinAlpha, 0, cosAlpha, 0,
    0, 0, 0, 1;

    return res;
}
void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    frameBuffer.setZero();
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        VertexAttributes attribute;
        attribute.position = uniform.cam * va.position;
        return attribute;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: build the vertex attributes from vertices and facets
    
    Vector3d x;
    Vector3d y;
    Vector3d z;
    Vector3i index;
    
    for (int i = 0; i<facets.rows(); ++i) {
        index << facets.row(i).transpose();
        x << vertices.row(index(0)).transpose();
        y << vertices.row(index(1)).transpose();
        z << vertices.row(index(2)).transpose();
        const VertexAttributes xa(x[0],x[1],x[2]);
        const VertexAttributes ya(y[0],y[1],y[2]);
        const VertexAttributes za(z[0],z[1],z[2]);

        vertex_attributes.push_back(xa);
        vertex_attributes.push_back(ya);
        vertex_attributes.push_back(za);
    }
    

    
    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}



void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    frameBuffer.setZero();
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    Matrix4d trafo = compute_rotation(alpha);

    

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform)
    {
        VertexAttributes attribute;
        attribute.position = uniform.cam * va.position;
        return attribute;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform)
    {
        // TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous)
    {
        // TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;
    Vector3d x;
    Vector3d y;
    Vector3d z;
    Vector3i index;
    for (int i = 0; i<facets.rows(); ++i) {
        index << facets.row(i).transpose();
        x << vertices.row(index(0)).transpose();
        y << vertices.row(index(1)).transpose();
        z << vertices.row(index(2)).transpose();
        const VertexAttributes xa(x[0],x[1],x[2]);
        const VertexAttributes ya(y[0],y[1],y[2]);
        const VertexAttributes za(z[0],z[1],z[2]);

        vertex_attributes.push_back(xa);
        vertex_attributes.push_back(ya);
        vertex_attributes.push_back(ya);
        vertex_attributes.push_back(za);
        vertex_attributes.push_back(za);
        vertex_attributes.push_back(xa);
    }
    
    
    uniform.cam *= trafo;
    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}


void get_shading_program(Program &program)
{
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: transform the position and the normal
        //TODO: compute the correct lighting
        return va;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: create the correct fragment
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: implement the depth check
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);
    Eigen::Matrix4d trafo = compute_rotation(alpha);

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: compute the normals
    //TODO: set material colors

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);

    //TODO: compute the vertex normals as vertex normal average

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: create vertex attributes
    //TODO: set material colors

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    
    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    //TODO: add the animation
    
    GifWriter g;
    GifBegin(&g, "wireframe.gif", frameBuffer.rows(), frameBuffer.cols(), 20);

    for (float i = 0; i < 20; i += EIGEN_PI / 10)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        wireframe_render(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), 20);
    }

    GifEnd(&g);

    return 0;
}
