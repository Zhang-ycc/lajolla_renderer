Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector = normalize(dir_in + dir_out);

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real n_in = fabs(dot(frame.n, dir_in));
    Real n_out = fabs(dot(frame.n, dir_out));
    Real h_out = fabs(dot(half_vector, dir_out));

    Real f_d90 = 1./2 + 2 * roughness * pow(h_out, 2);
    Real f_dIn = 1 + (f_d90 - 1) * pow(1 - n_in, 5);
    Real f_dOut = 1 + (f_d90 - 1) * pow(1 - n_out, 5);

    Real f_ss90 = roughness * pow(h_out, 2);
    Real f_ssIn = 1 + (f_ss90 - 1) * pow(1 - n_in, 5);
    Real f_ssOut = 1 + (f_ss90 - 1) * pow(1 - n_out, 5);

    Spectrum f_baseDiffuse = base_color / c_PI * f_dIn * f_dOut * n_out;
    Spectrum f_subsurface = 1.25 * base_color / c_PI * (f_ssIn * f_ssOut * (1 / (n_in + n_out) - 0.5) + 0.5) * n_out;
    Spectrum f_diffuse = (1 - subsurface) * f_baseDiffuse + subsurface * f_subsurface;

    return f_diffuse;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    // return Real(0);
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */};
    // return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
