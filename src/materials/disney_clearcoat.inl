#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_out = fabs(dot(half_vector, dir_out));

    Real n = 1.5;
    Spectrum R_0 = make_const_spectrum(pow(n - 1, 2) / pow(n + 1, 2));
    Spectrum F_c = R_0 + (1 - R_0) * pow(1 - h_out, 5);

    Real a_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Vector3 h_l = to_local(frame, half_vector);
    Real D_c = (a_g * a_g - 1) / (c_PI * log(a_g * a_g) * (1 + (a_g * a_g - 1) * h_l.z * h_l.z));

    Vector3 w_in_l = to_local(frame, dir_in);
    Vector3 w_out_l = to_local(frame, dir_out);
    Real V_in = (sqrt(1 + (pow(w_in_l.x * 0.25, 2) + pow(w_in_l.y * 0.25, 2)) / pow(w_in_l.z, 2)) - 1) / 2;
    Real V_out = (sqrt(1 + (pow(w_out_l.x * 0.25, 2) + pow(w_out_l.y * 0.25, 2)) / pow(w_out_l.z, 2)) - 1) / 2;
    Real G_in = 1 / (1 + V_in);
    Real G_out = 1 / (1 + V_out);
    Real G_c = G_in * G_out;

    Spectrum f_clearcoat = F_c * D_c * G_c / (4 * fabs(dot(frame.n, dir_in)));

    return f_clearcoat;

    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_out = fabs(dot(half_vector, dir_out));

    Real a_g = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Vector3 h_l = to_local(frame, half_vector);
    Real D_c = (a_g * a_g - 1) / (c_PI * log(a_g * a_g) * (1 + (a_g * a_g - 1) * h_l.z * h_l.z));

    return D_c * fabs(dot(frame.n, half_vector)) / (4 * h_out);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real a = (1 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real h_e = acos(sqrt((1 - pow(a * a, 1 - rnd_param_uv.x)) / (1 - a * a)));
    Real h_a = 2 * c_PI * rnd_param_uv.y;
    Vector3 h_l = Vector3(sin(h_e) * cos(h_a), sin(h_e) * sin(h_a), cos(h_e));
    Vector3 half_vector = to_world(frame, h_l);
    Vector3 reflected = normalize(- dir_in + 2 * dot(dir_in, half_vector) * half_vector);

    return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, Real(1) /* roughness */};

    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}