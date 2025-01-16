#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    //Glass
    DisneyGlass disneyGlass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    Spectrum f_glass = this->operator()(disneyGlass);

    if (dot(dir_in, vertex.geometric_normal) <= 0 || !reflect) {
        return (1 - metallic) * specularTransmission * f_glass;
    }

    DisneyDiffuse disneyDiffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    Spectrum f_diffuse = this->operator()(disneyDiffuse);

    DisneyClearcoat disneyClearcoat = {bsdf.clearcoat_gloss};
    Spectrum f_clearcoat = this->operator()(disneyClearcoat);

    DisneySheen disneySheen = {bsdf.base_color, bsdf.sheen_tint};
    Spectrum f_sheen = this->operator()(disneySheen);
    
    //Metal
    Spectrum f_metal;
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        f_metal = make_zero_spectrum();
    } else {
        Spectrum white = make_const_spectrum(Real(1.));
        Spectrum C_tint = luminance(base_color) > 0 ? base_color / luminance(base_color) : white;
        Spectrum K_s = (1 - specular_tint) + specular_tint * C_tint;
        Spectrum R_0 = make_const_spectrum(pow(eta - 1, 2) / pow(eta + 1, 2));
        Spectrum C_0 = specular * R_0 * (1 - metallic) * K_s + metallic * base_color;
        Vector3 half_vector = normalize(dir_in + dir_out);
        Spectrum F_m = C_0 + (1 - C_0) * pow(1 - dot(half_vector, dir_out), 5);

        roughness = std::clamp(roughness, Real(0.01), Real(1));
        Real a_min = 0.0001;
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

        f_metal = F_m * D_m * G_m / (4 * fabs(dot(frame.n, dir_in)));
    }


    return (1 - specularTransmission) * (1 - metallic) * f_diffuse
        + (1 - metallic) * sheen * f_sheen
        + (1 - specularTransmission * (1 - metallic)) * f_metal
        + 0.25 * clearcoat * f_clearcoat
        + (1 - metallic) * specularTransmission * f_glass;
    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    DisneyGlass disneyGlass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};

    if (dot(dir_in, vertex.geometric_normal) <= 0 || !reflect) {
        return this->operator()(disneyGlass);
    }

    DisneyDiffuse disneyDiffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyClearcoat disneyClearcoat = {bsdf.clearcoat_gloss};
    DisneyMetal disneyMetal = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic};

    Real diffuseWeight = (1 - metallic) * (1 - specularTransmission);
    Real metalWeight = (1 - specularTransmission * (1 - metallic));
    Real glassWeight = (1 - metallic) * specularTransmission;
    Real clearcoatWeight = 0.25 * clearcoat;
    Real totalWeight = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;

    return diffuseWeight / totalWeight * this->operator()(disneyDiffuse)
        + metalWeight / totalWeight * this->operator()(disneyMetal)
        + glassWeight / totalWeight * this->operator()(disneyGlass)
        + clearcoat / totalWeight * this->operator()(disneyClearcoat);
    // return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    DisneyGlass disneyGlass = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    if (dot(dir_in, vertex.geometric_normal) <= 0) {
        return this->operator()(disneyGlass);
    }

    Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real diffuseWeight = (1 - metallic) * (1 - specularTransmission);
    Real metalWeight = (1 - specularTransmission * (1 - metallic));
    Real glassWeight = (1 - metallic) * specularTransmission;
    Real clearcoatWeight = 0.25 * clearcoat;
    Real totalWeight = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;
    diffuseWeight /= totalWeight;
    metalWeight /= totalWeight;
    glassWeight /= totalWeight;
    clearcoatWeight /= totalWeight;

    Real w = rnd_param_w;
    if (w < diffuseWeight) {
        DisneyDiffuse disneyDiffuse = {bsdf.base_color, bsdf.roughness, bsdf.subsurface};
        return this->operator()(disneyDiffuse);
    } else if (w < diffuseWeight + metalWeight) {
        DisneyMetal disneyMetal = {bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
        return this->operator()(disneyMetal);
    } else if (w < diffuseWeight + metalWeight + glassWeight) {
        Real rescaleW = w / (diffuseWeight + metalWeight + glassWeight);
        return sample_bsdf(disneyGlass, dir_in, vertex, texture_pool, rnd_param_uv, rescaleW, dir);
    } else if (w < diffuseWeight + metalWeight + glassWeight + clearcoatWeight) {
        DisneyClearcoat disneyClearcoat = {bsdf.clearcoat_gloss};
        return this->operator()(disneyClearcoat);
    }

    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
