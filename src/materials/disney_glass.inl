#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw

    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = fmax(a_min, roughness * roughness / aspect);
    Real a_y = fmax(a_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);
    Real D = 1 / (c_PI * a_x * a_y * pow(pow(h_l.x / a_x, 2) + pow(h_l.y / a_y, 2) + h_l.z * h_l.z, 2));

    Vector3 w_in_l = to_local(frame, dir_in);
    Vector3 w_out_l = to_local(frame, dir_out);
    Real V_in = (sqrt(1 + (pow(w_in_l.x * a_x, 2) + pow(w_in_l.y * a_y, 2)) / pow(w_in_l.z, 2)) - 1) / 2;
    Real V_out = (sqrt(1 + (pow(w_out_l.x * a_x, 2) + pow(w_out_l.y * a_y, 2)) / pow(w_out_l.z, 2)) - 1) / 2;
    Real G_in = 1 / (1 + V_in);
    Real G_out = 1 / (1 + V_out);
    Real G = G_in * G_out;

    Real F = fresnel_dielectric(h_dot_in, eta);
    if (reflect) {
        return base_color * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        
        return sqrt(base_color) * ((1 - F) * D * G * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
    }

    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = fmax(a_min, roughness * roughness / aspect);
    Real a_y = fmax(a_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);
    Real D = 1 / (c_PI * a_x * a_y * pow(pow(h_l.x / a_x, 2) + pow(h_l.y / a_y, 2) + h_l.z * h_l.z, 2));

    Vector3 w_in_l = to_local(frame, dir_in);
    Real V_in = (sqrt(1 + (pow(w_in_l.x * a_x, 2) + pow(w_in_l.y * a_y, 2)) / pow(w_in_l.z, 2)) - 1) / 2;
    Real G_in = 1 / (1 + V_in);

    
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }

    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a_min = 0.0001;
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Vector3 local_dir_in = to_local(frame, dir_in);

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real a_x = fmax(a_min, roughness * roughness / aspect);
    Real a_y = fmax(a_min, roughness * roughness * aspect);
    Vector3 local_micro_normal = sample_visible_normals(local_dir_in, a_x, a_y, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }

    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
