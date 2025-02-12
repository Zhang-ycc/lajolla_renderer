#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    const Medium &media = scene.media[scene.camera.medium_id];

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (vertex_) {
        PathVertex vertex = *vertex_;
        Spectrum sigma_a = get_sigma_a(media, vertex.position);
        Real t = distance(ray.org, vertex.position);
        Spectrum transmittance = exp(-sigma_a * t);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return transmittance * Le;
    }

    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this、、
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    const Medium &media = scene.media[scene.camera.medium_id];
    Spectrum sigma_a = get_sigma_a(media, ray.org);
    Spectrum sigma_s = get_sigma_s(media, ray.org);
    Real sigma_t = (sigma_a + sigma_s)[0];

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t;

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    Real t_hit = vertex_ ? distance(ray.org, vertex_->position) : infinity<Real>();

    if (t < t_hit) {
        Real trans_pdf = exp(-sigma_t * t) * sigma_t;
        Real transmittance = exp(-sigma_t * t);
        
        // compute L_s1 using Monte Carlo sampling
        Spectrum p = ray.org + t * ray.dir;

        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);
        Spectrum dir_light = normalize(point_on_light.position - p);

        Spectrum rho = eval(get_phase_function(media), -ray.dir, dir_light);
        Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real tr = exp(- sigma_t * distance(point_on_light.position, p));
        
        Real G = 0;
        Ray shadow_ray{p, dir_light, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(point_on_light.position, p)};
        if (!occluded(scene, shadow_ray)) {
            G = max(-dot(dir_light, point_on_light.normal), Real(0)) / distance_squared(point_on_light.position, p);
        }

        Spectrum L_s1_estimate = rho * Le * tr * G;
        Real L_s1_pdf = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, p, scene);

        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf);
    } else {
        PathVertex vertex = *vertex_;
        Real trans_pdf = exp(-sigma_t * t_hit);
        Real transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return transmittance / trans_pdf * Le;
    }
}

int update_medium(PathVertex isect, Ray ray, int medium_id) {
    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // At medium transition. Update medium.
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            medium_id = isect.exterior_medium_id;
        } else {
            medium_id = isect.interior_medium_id;
        }
    }
    return medium_id;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;

    while (true) {
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_ ? distance(ray.org, vertex_->position) : infinity<Real>();

        const Medium &media = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(media, ray.org);
        Spectrum sigma_s = get_sigma_s(media, ray.org);
        Real sigma_t = (sigma_a + sigma_s)[0];
        
        Real transmittance = 1.0;
        Real trans_pdf = 1.0;
        if (current_medium_id != -1) {
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            // compute transmittance and trans_pdf
            // if t < t_hit, set scatter = True
            if (t < t_hit) {
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                ray.org = vertex_->position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
            }
        }

        current_path_throughput *= (transmittance / trans_pdf);
        
        if (!scatter && vertex_) {
            // reach a surface, include emission
            PathVertex vertex = *vertex_;
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[vertex.shape_id])) {
                Le = emission(vertex, -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }

        int max_depth = scene.options.max_depth;
        if (bounces == max_depth - 1 && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                ray.org = vertex.position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phaseFunction = get_phase_function(media);
            std::optional<Spectrum> next_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_param);
            if (!next_dir) {
                // phase function sampling failed. Abort the loop.
                break;
            }

            Spectrum rho = eval(phaseFunction, -ray.dir, *next_dir);
            current_path_throughput *= (rho / pdf_sample_phase(phaseFunction, -ray.dir, *next_dir)) * sigma_s;
            // update ray.dir
            ray.dir = *next_dir;
        } else {
            // Hit a surface -- don’t need to deal with this yet
            break;
        }

        Real rr_prob = 1.0;
        int rr_depth = scene.options.rr_depth;
        if (bounces >= rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        
        bounces += 1;
    }

    return radiance;
}

Spectrum next_event_estimation_volume(const Scene &scene, Spectrum p, int current_medium_id, Ray ray, int bounces, pcg32_state &rng) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Spectrum dir_light = normalize(p_prime.position - p);
    
    // // Compute transmittance to light. Skip through index-matching shapes.
    Real T_light = 1.0;
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.0; // for multiple importance sampling
    int max_depth = scene.options.max_depth;
    Spectrum ori_p = p;

    while (true) {
        Ray shadow_ray{p, dir_light,  get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        std::optional<PathVertex> isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }
        const Medium &media = scene.media[shadow_medium_id];
        Spectrum sigma_a = get_sigma_a(media, p);
        Spectrum sigma_s = get_sigma_s(media, p);
        Real sigma_t = (sigma_a + sigma_s)[0];
        
        // Account for the transmittance to next_t
        if (shadow_medium_id != -1) {
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // Something is blocking: is it an opaque surface?
            if (isect->material_id >= 0) {
                // we’re blocked
                return make_zero_spectrum();
            }
            // otherwise, it’s an index-matching surface and we want to pass through 
            //-- this introduces one extra connection vertex
            shadow_bounces += 1;
            if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(*isect, shadow_ray, shadow_medium_id);
            p = p + next_t * dir_light;
        }
    }

    if (T_light > 0) {
        // Compute T_light * G * rho * L & pdf_nee
        PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);

        Real G = max(-dot(dir_light, p_prime.normal), Real(0)) / distance_squared(p_prime.position, ori_p);
        Spectrum rho = eval(phaseFunction, -ray.dir, dir_light);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, ori_p, scene);
        Spectrum contrib = T_light * G * rho * L / pdf_nee;

        Real pdf_phase = pdf_sample_phase(phaseFunction, -ray.dir, dir_light) * G * p_trans_dir;
        // power heuristics
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

        return w * contrib;
    }

    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0; //in solid angle measure
    Spectrum nee_p_cache;
    Real multi_trans_pdf = 1;
    Real never_scatter = true;

    while (true) {
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_ ? distance(ray.org, vertex_->position) : infinity<Real>();

        const Medium &media = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(media, ray.org);
        Spectrum sigma_s = get_sigma_s(media, ray.org);
        Real sigma_t = (sigma_a + sigma_s)[0];
        
        Real transmittance = 1.0;
        Real trans_pdf = 1.0;
        if (current_medium_id != -1) {
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            // compute transmittance and trans_pdf
            // if t < t_hit, set scatter = True
            if (t < t_hit) {
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                ray.org = vertex_->position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
            }
            multi_trans_pdf *= trans_pdf;
        }

        current_path_throughput *= (transmittance / trans_pdf);
        
        // If we reach a surface and didn’t scatter, include the emission.
        if (!scatter && vertex_) {
            if (never_scatter) {
                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                PathVertex vertex = *vertex_;
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                PathVertex vertex = *vertex_;
                if (is_light(scene.shapes[vertex.shape_id])) {
                    // Need to account for next event estimation
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    assert(light_id >= 0);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    // Note that pdf_nee needs to account for the path vertex that issued
                    // next event estimation potentially many bounces ago.
                    // The vertex position is stored in nee_p_cache.
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene);
                    // The PDF for sampling the light source using phase function sampling + transmittance sampling
                    // The directional sampling pdf was cached in dir_pdf in solid angle measure.
                    // The transmittance sampling pdf was cached in multi_trans_pdf.
                    Real G = max(-dot(ray.dir, vertex.geometric_normal), Real(0)) / distance_squared(vertex.position, nee_p_cache);
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    // current_path_throughput already accounts for transmittance.
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        int max_depth = scene.options.max_depth;
        if (bounces == max_depth - 1 && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                ray.org = vertex.position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            Spectrum nee_contrib = next_event_estimation_volume(scene, ray.org, current_medium_id, ray, bounces, rng);
            radiance += current_path_throughput * nee_contrib * sigma_s;

            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phaseFunction = get_phase_function(media);
            std::optional<Spectrum> next_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_param);
            if (!next_dir) {
                // phase function sampling failed. Abort the loop.
                break;
            }

            Spectrum rho = eval(phaseFunction, -ray.dir, *next_dir);
            Real pdf_phase_sample = pdf_sample_phase(phaseFunction, -ray.dir, *next_dir);
            current_path_throughput *= (rho / pdf_phase_sample) * sigma_s;
            // update ray.dir
            ray.dir = *next_dir;

            // store pdf, nee_p_cache
            dir_pdf = pdf_phase_sample;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;
        } else {
            // Hit a surface -- don’t need to deal with this yet
            break;
        }

        Real rr_prob = 1.0;
        int rr_depth = scene.options.rr_depth;
        if (bounces >= rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        
        bounces += 1;
        never_scatter = false;
    }

    return radiance;
}

Spectrum next_event_estimation_surface(const Scene &scene, Spectrum p, int current_medium_id, Ray ray, int bounces, PathVertex vertex, pcg32_state &rng) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Spectrum dir_light = normalize(p_prime.position - p);
    
    // Compute transmittance to light. Skip through index-matching shapes.
    Real T_light = 1.0;
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.0; // for multiple importance sampling
    int max_depth = scene.options.max_depth;
    Spectrum ori_p = p;

    while (true) {
        Ray shadow_ray{p, dir_light,  get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        std::optional<PathVertex> isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }
        const Medium &media = scene.media[shadow_medium_id];
        Spectrum sigma_a = get_sigma_a(media, p);
        Spectrum sigma_s = get_sigma_s(media, p);
        Real sigma_t = (sigma_a + sigma_s)[0];
        
        // Account for the transmittance to next_t
        if (shadow_medium_id != -1) {
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // Something is blocking: is it an opaque surface?
            if (isect->material_id >= 0) {
                // we’re blocked
                return make_zero_spectrum();
            }
            // otherwise, it’s an index-matching surface and we want to pass through 
            //-- this introduces one extra connection vertex
            shadow_bounces += 1;
            if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(*isect, shadow_ray, shadow_medium_id);
            p = p + next_t * dir_light;
        }
    }

    if (T_light > 0) {
        // surface
        const Material &mat = scene.materials[vertex.material_id];
        Spectrum f = eval(mat, -ray.dir, dir_light, vertex, scene.texture_pool);

        Real G = max(-dot(dir_light, p_prime.normal), Real(0)) / distance_squared(p_prime.position, ori_p);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, ori_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;

        Real pdf_bsdf_sample = pdf_sample_bsdf(mat, -ray.dir, dir_light, vertex, scene.texture_pool) * G * p_trans_dir;
        // power heuristics
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf_sample * pdf_bsdf_sample);

        return w * contrib;
    }

    return make_zero_spectrum();
}


// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0; //in solid angle measure
    Spectrum nee_p_cache;
    Real multi_trans_pdf = 1.0;
    Real never_scatter = true;
    Real eta_scale = 1.0;

    while (true) {
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_ ? distance(ray.org, vertex_->position) : infinity<Real>();

        const Medium &media = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(media, ray.org);
        Spectrum sigma_s = get_sigma_s(media, ray.org);
        Real sigma_t = (sigma_a + sigma_s)[0];
        
        Real transmittance = 1.0;
        Real trans_pdf = 1.0;
        if (current_medium_id != -1) {
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            // compute transmittance and trans_pdf
            // if t < t_hit, set scatter = True
            if (t < t_hit) {
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                scatter = true;
                ray.org = ray.org + t * ray.dir;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                ray.org = vertex_->position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
            }
            multi_trans_pdf *= trans_pdf;
        }

        current_path_throughput *= (transmittance / trans_pdf);
        
        // If we reach a surface and didn’t scatter, include the emission.
        if (!scatter && vertex_) {
            ray.org = vertex_->position;
            ray.tnear= get_intersection_epsilon(scene);
            
            if (never_scatter) {
                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                PathVertex vertex = *vertex_;
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                PathVertex vertex = *vertex_;
                if (is_light(scene.shapes[vertex.shape_id])) {
                    // Need to account for next event estimation
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    assert(light_id >= 0);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    // Note that pdf_nee needs to account for the path vertex that issued
                    // next event estimation potentially many bounces ago.
                    // The vertex position is stored in nee_p_cache.
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene);
                    // The PDF for sampling the light source using phase function sampling + transmittance sampling
                    // The directional sampling pdf was cached in dir_pdf in solid angle measure.
                    // The transmittance sampling pdf was cached in multi_trans_pdf.
                    Real G = max(-dot(ray.dir, vertex.geometric_normal), Real(0)) / distance_squared(vertex.position, nee_p_cache);
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    // current_path_throughput already accounts for transmittance.
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        int max_depth = scene.options.max_depth;
        if (bounces == max_depth - 1 && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                ray.org = vertex.position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            Spectrum nee_contrib = next_event_estimation_volume(scene, ray.org, current_medium_id, ray, bounces, rng);
            radiance += current_path_throughput * nee_contrib * sigma_s;

            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phaseFunction = get_phase_function(media);
            std::optional<Spectrum> next_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_param);
            if (!next_dir) {
                // phase function sampling failed. Abort the loop.
                break;
            }

            Spectrum rho = eval(phaseFunction, -ray.dir, *next_dir);
            Real pdf_phase_sample = pdf_sample_phase(phaseFunction, -ray.dir, *next_dir);

            current_path_throughput *= (rho / pdf_phase_sample) * sigma_s;
            // update ray.dir
            ray.dir = *next_dir;

            // store pdf, nee_p_cache
            dir_pdf = pdf_phase_sample;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;
        } else if (vertex_) {
            // Hit a surface -- don’t need to deal with this yet
            PathVertex vertex = *vertex_;
            Spectrum nee_contrib = next_event_estimation_surface(scene, ray.org, current_medium_id, ray, bounces, vertex,rng);
            radiance += current_path_throughput * nee_contrib;
            
            const Material &mat = scene.materials[vertex.material_id];

            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, -ray.dir, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;

            // Update ray differentials & eta_scale
            if (bsdf_sample.eta != 0) {
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
            }

            Spectrum f = eval(mat, -ray.dir, dir_bsdf, vertex, scene.texture_pool);
            Real pdf_bsdf_sample = pdf_sample_bsdf(mat, -ray.dir, dir_bsdf, vertex, scene.texture_pool);

            current_path_throughput *= (f / pdf_bsdf_sample);
            // update ray.dir
            // Trace a ray towards bsdf_dir. Note that again we have
            // to have an "epsilon" tnear to prevent self intersection.
            ray = Ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};

            current_medium_id = update_medium(vertex, ray, current_medium_id);
            // store pdf, nee_p_cache
            dir_pdf = pdf_bsdf_sample;
            nee_p_cache = ray.org;
            multi_trans_pdf = 1;

            // break;
        }

        Real rr_prob = 1.0;
        int rr_depth = scene.options.rr_depth;
        if (bounces >= rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        
        bounces += 1;
        never_scatter = false;
    }

    return radiance;
}

Spectrum next_event_estimation(const Scene &scene, Spectrum p, int current_medium_id, Ray ray, int bounces, pcg32_state &rng, bool is_surface, PathVertex vertex = {}) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Spectrum dir_light = normalize(p_prime.position - p);
    
    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light = make_const_spectrum(1);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Spectrum p_trans_dir = make_const_spectrum(1); // for multiple importance sampling  
    Spectrum p_trans_nee = make_const_spectrum(1);
    int max_depth = scene.options.max_depth;
    int max_null_collisions = scene.options.max_null_collisions;
    Spectrum ori_p = p;

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    while (true) {
        const Medium &media = scene.media[shadow_medium_id];

        Ray shadow_ray{p, dir_light, get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            if (isect->material_id >= 0) {
                // we’re blocked
                return make_zero_spectrum();
            }
            next_t = distance(p, isect->position);
        }

        // Account for the transmittance to next_t(including p_prime!)
        if (shadow_medium_id != -1) {
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;
            Ray iter_ray{p, dir_light};
            Spectrum majorant = get_majorant(media, shadow_ray);

            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= max_null_collisions) {
                    break;
                }

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;
                // Update accumulated distance
                accum_t = min(accum_t + t, next_t);
                if (t < dt){
                    iter_ray.org += t * iter_ray.dir;

                    Spectrum sigma_a = get_sigma_a(media, iter_ray.org);
                    Spectrum sigma_s = get_sigma_s(media, iter_ray.org);
                    Spectrum sigma_t = sigma_a + sigma_s;
                    Spectrum sigma_n = majorant - sigma_t;

                    // didn’t hit the surface, so this is a null-scattering event
                    T_light *= exp(-majorant * t) * sigma_n / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        // optimization for places where sigma_n = 0
                        break;
                    }
                } else {
                    // hit the surface
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    // iter_ray.org = isect->position;
                    break;
                }
                iteration += 1;
            }
        }
        

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // Something is blocking: is it an opaque surface?
            if (isect->material_id >= 0) {
                // we’re blocked
                return make_zero_spectrum();
            }
            // otherwise, it’s an index-matching surface and we want to pass through 
            //-- this introduces one extra connection vertex
            shadow_bounces += 1;
            if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(*isect, shadow_ray, shadow_medium_id);
            p = p + next_t * dir_light;
        }
    }

    if (max(T_light) > 0) {
        if (!is_surface) {
            // Compute T_light * G * rho * L & pdf_nee
            PhaseFunction phaseFunction = get_phase_function(scene.media[current_medium_id]);

            Real G = max(-dot(dir_light, p_prime.normal), Real(0)) / distance_squared(p_prime.position, ori_p);
            Spectrum rho = eval(phaseFunction, -ray.dir, dir_light);
            Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
            Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, ori_p, scene) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * rho * L / pdf_nee;

            Real pdf_phase = pdf_sample_phase(phaseFunction, -ray.dir, dir_light) * G * avg(p_trans_dir);
            // power heuristics
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

            return w * contrib;
        } else {
            // surface
            const Material &mat = scene.materials[vertex.material_id];
            Spectrum f = eval(mat, -ray.dir, dir_light, vertex, scene.texture_pool);

            Real G = max(-dot(dir_light, p_prime.normal), Real(0)) / distance_squared(p_prime.position, ori_p);
            Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
            Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, ori_p, scene) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;

            Real pdf_bsdf_sample = pdf_sample_bsdf(mat, -ray.dir, dir_light, vertex, scene.texture_pool) * G * avg(p_trans_dir);
            // power heuristics
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf_sample * pdf_bsdf_sample);

            return w * contrib;
        }
    }

    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
        int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    
    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0; //in solid angle measure
    Spectrum nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_nee_pdf = make_const_spectrum(1);
    Real never_scatter = true;
    Real eta_scale = 1.0;

    while (true) {
        bool scatter = false;

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_ ? distance(ray.org, vertex_->position) : infinity<Real>();

        const Medium &media = scene.media[current_medium_id];
        
        Spectrum transmittance = make_const_spectrum(1);
        // Spectrum trans_pdf = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1); // PDF for free-flight sampling
        Spectrum trans_nee_pdf = make_const_spectrum(1); // PDF for next event estimation

        int max_null_collisions = scene.options.max_null_collisions;

        if (current_medium_id != -1) {
            // Sample a channel for sampling
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;
            Spectrum majorant = get_majorant(media, ray);

            while (true) {

                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= max_null_collisions) {
                    break;
                }

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                // Update accumulated distance
                accum_t = min(accum_t + t, t_hit);
                if (t < dt){
                    ray.org = ray.org + t * ray.dir;

                    Spectrum sigma_a = get_sigma_a(media, ray.org);
                    Spectrum sigma_s = get_sigma_s(media, ray.org);
                    Spectrum sigma_t = sigma_a + sigma_s;
                    Spectrum sigma_n = majorant - sigma_t;

                    // haven’t reached the surface
                    // sample from real/fake particle events
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        // hit a "real" particle
                        scatter = true;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        // don’t need to account for trans_nee_pdf since we scatter
                        ray.tnear= get_intersection_epsilon(scene);
                        break;
                    } else {
                        // hit a "fake" particle
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    // reach the surface
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    ray.org = vertex_->position;
                    // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                    ray.tnear= get_intersection_epsilon(scene);
                    break;
                }
                iteration += 1;
            }
            multi_trans_pdf *= trans_dir_pdf;
            multi_nee_pdf *= trans_nee_pdf;
        }

        current_path_throughput *= (transmittance / avg(trans_dir_pdf));
        
        // If we reach a surface and didn’t scatter, include the emission.
        if (!scatter && vertex_) {
            ray.org = vertex_->position;
            ray.tnear= get_intersection_epsilon(scene);
            
            if (never_scatter) {
                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                PathVertex vertex = *vertex_;
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[vertex.shape_id])) {
                    Le = emission(vertex, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le;
            } else {
                PathVertex vertex = *vertex_;
                if (is_light(scene.shapes[vertex.shape_id])) {
                    // Need to account for next event estimation
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    assert(light_id >= 0);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    // Note that pdf_nee needs to account for the path vertex that issued
                    // next event estimation potentially many bounces ago.
                    // The vertex position is stored in nee_p_cache.
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene) * avg(multi_nee_pdf);
                    // The PDF for sampling the light source using phase function sampling + transmittance sampling
                    // The directional sampling pdf was cached in dir_pdf in solid angle measure.
                    // The transmittance sampling pdf was cached in multi_trans_pdf.
                    Real G = max(-dot(ray.dir, vertex.geometric_normal), Real(0)) / distance_squared(vertex.position, nee_p_cache);
                    Real dir_pdf_ = dir_pdf * avg(multi_trans_pdf) * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    // current_path_throughput already accounts for transmittance.
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
                }
            }
        }

        int max_depth = scene.options.max_depth;
        if (bounces == max_depth - 1 && max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                // index-matching interface, skip through it
                ray.org = vertex.position;
                // ray origin from surface, have an "epsilon" tnear to prevent self intersection.
                ray.tnear= get_intersection_epsilon(scene);
                current_medium_id = update_medium(vertex, ray, current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if (scatter) {
            Spectrum sigma_s = get_sigma_s(media, ray.org);
            Spectrum nee_contrib = next_event_estimation(scene, ray.org, current_medium_id, ray, bounces, rng, false);
            // Spectrum nee_contrib = next_event_estimation_volume(scene, ray.org, current_medium_id, ray, bounces, rng);
            radiance += current_path_throughput * nee_contrib * sigma_s;

            Vector2 rnd_param{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phaseFunction = get_phase_function(media);
            std::optional<Spectrum> next_dir = sample_phase_function(phaseFunction, -ray.dir, rnd_param);
            if (!next_dir) {
                // phase function sampling failed. Abort the loop.
                break;
            }

            Spectrum rho = eval(phaseFunction, -ray.dir, *next_dir);
            Real pdf_phase_sample = pdf_sample_phase(phaseFunction, -ray.dir, *next_dir);

            current_path_throughput *= (rho / pdf_phase_sample) * sigma_s;
            // update ray.dir
            ray.dir = *next_dir;

            // store pdf, nee_p_cache
            dir_pdf = pdf_phase_sample;
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1.0);
            multi_nee_pdf = make_const_spectrum(1.0);
        } else if (vertex_) {
            // Hit a surface -- don’t need to deal with this yet
            PathVertex vertex = *vertex_;
            Spectrum nee_contrib = next_event_estimation(scene, ray.org, current_medium_id, ray, bounces, rng, true, vertex);
            radiance += current_path_throughput * nee_contrib;
            
            const Material &mat = scene.materials[vertex.material_id];

            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ = sample_bsdf(mat, -ray.dir, vertex, scene.texture_pool, bsdf_rnd_param_uv, bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
            Vector3 dir_bsdf = bsdf_sample.dir_out;

            // Update ray differentials & eta_scale
            if (bsdf_sample.eta != 0) {
                eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
            }

            Spectrum f = eval(mat, -ray.dir, dir_bsdf, vertex, scene.texture_pool);
            Real pdf_bsdf_sample = pdf_sample_bsdf(mat, -ray.dir, dir_bsdf, vertex, scene.texture_pool);

            current_path_throughput *= (f / pdf_bsdf_sample);
            // update ray.dir
            // Trace a ray towards bsdf_dir. Note that again we have
            // to have an "epsilon" tnear to prevent self intersection.
            ray = Ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};

            current_medium_id = update_medium(vertex, ray, current_medium_id);
            // store pdf, nee_p_cache
            dir_pdf = pdf_bsdf_sample;
            nee_p_cache = ray.org;
            multi_trans_pdf = make_const_spectrum(1.0);
            multi_nee_pdf = make_const_spectrum(1.0);

            // break;
        }

        Real rr_prob = 1.0;
        int rr_depth = scene.options.rr_depth;
        if (bounces >= rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        
        bounces += 1;
        never_scatter = false;
    }

    return radiance;
}