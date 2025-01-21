#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;

    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_out = fabs(dot(half_vector, dir_out));

    Spectrum F_m = base_color + (1 - base_color) * pow(1 - h_out, 5);

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = fmax(a_min, roughness * roughness / aspect);
    Real a_y = fmax(a_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);

    Real D_m = 1 / (c_PI * a_x * a_y * pow(pow(h_l.x / a_x, 2) + pow(h_l.y / a_y, 2) + h_l.z * h_l.z, 2));

    Vector3 w_in_l = to_local(frame, dir_in);
    Vector3 w_out_l = to_local(frame, dir_out);
    Real V_in = (sqrt(1 + (pow(w_in_l.x * a_x, 2) + pow(w_in_l.y * a_y, 2)) / pow(w_in_l.z, 2)) - 1) / 2;
    Real V_out = (sqrt(1 + (pow(w_out_l.x * a_x, 2) + pow(w_out_l.y * a_y, 2)) / pow(w_out_l.z, 2)) - 1) / 2;
    Real G_in = 1 / (1 + V_in);
    Real G_out = 1 / (1 + V_out);

    Real G_m = G_in * G_out;

    Spectrum f_metal = F_m * D_m * G_m / (4 * fabs(dot(frame.n, dir_in)));

    return f_metal;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;

    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 half_vector = normalize(dir_in + dir_out);

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = std::max(a_min, roughness * roughness / aspect);
    Real a_y = std::max(a_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);

    Real D_m = 1 / (c_PI * a_x * a_y * pow(h_l.x * h_l.x / a_x / a_x + h_l.y * h_l.y / a_y / a_y + h_l.z * h_l.z, 2));

    Vector3 w_in_l = to_local(frame, dir_in);
    Real V_in = (sqrt(1 + (pow(w_in_l.x * a_x, 2) + pow(w_in_l.y * a_y, 2)) / pow(w_in_l.z, 2)) - 1) / 2;
    Real G_in = 1 / (1 + V_in);

    return (D_m * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    // Vector3 local_dir_in = to_local(frame, dir_in);
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;
    
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = fmax(a_min, roughness * roughness / aspect);
    Real a_y = fmax(a_min, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals(local_dir_in, a_x, a_y, rnd_param_uv);
        
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };

    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
