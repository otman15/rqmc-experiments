// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <gtest/gtest.h>

#include <pbrt/util/image.h>
#include <pbrt/util/noise.h>
#include <pbrt/util/rng.h>
#include <pbrt/util/progressreporter.h>
#include <pbrt/util/vecmath.h>

#include <cmath>
#include <vector>

using namespace pbrt;

TEST(Noise, StochasticMatch) {
    RNG rng;
    for (int i = 0; i < 100; ++i) {
        Float sum = 0;
        int n = 1024*1024;
        auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
        Point3f p(u(), u(), u());
        Float ref = Noise(p);
        for (int j = 0; j < n; ++j)
            sum += NoiseStoch(p, rng.Uniform<Float>());
        EXPECT_LT(std::abs(ref - sum/n), .007);
    }
}

TEST(Noise, FBmStochMatch) {
    RNG rng;
    int maxOctaves = 8;
    Vector3f dpdx(1e-8, 0, 0), dpdy(0, 1e-8, 0);
    for (int i = 0; i < 10; ++i) {
        Float sum = 0;
        int n = 1024*1024;
        auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
        Point3f p(u(), u(), u());

        Float ref = FBm(p, dpdx, dpdy, 0.6, maxOctaves);
        for (int j = 0; j < n; ++j)
            sum += FBmStoch(p, dpdx, dpdy, 0.6, maxOctaves, rng.Uniform<Float>());

        EXPECT_LT(std::abs(ref - sum/n), .007);
    }
}

TEST(Noise, Perf) {
    RNG rng;
    std::vector<Point3f> p;
    auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
    for (int i = 0; i < 100; ++i)
        p.push_back(Point3f(u(), u(), u()));

    int nbench = 1000000;
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += Noise(p[j]);
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Normal: %.3f (sum %f)\n", elapsed, sum);
    }
    {
        float u = rng.Uniform<Float>();
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += NoiseStoch(p[j], u);
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Stochastic: %.3f (sum %f)\n", elapsed, sum);
    }
}

TEST(Noise, FbmStochPerf) {
    RNG rng;
    std::vector<Point3f> p;
    auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
    for (int i = 0; i < 100; ++i)
        p.push_back(Point3f(u(), u(), u()));

    int maxOctaves = 8;
    Vector3f dpdx(1e-8, 0, 0), dpdy(0, 1e-8, 0);
    Float omega = 0.6;
    int nbench = 100000;
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += FBm(p[j], dpdx, dpdy, omega, maxOctaves);
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Normal: %.3f (sum %f)\n", elapsed, sum);
    }
    {
        float u = rng.Uniform<Float>();
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += FBmStoch(p[j], dpdx, dpdy, omega, maxOctaves, rng.Uniform<Float>());
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Stochastic: %.3f (sum %f)\n", elapsed, sum);
    }
}

TEST(Noise, LewisImage) {
    Image image(PixelFormat::Float, {512, 512}, {"Y"});

    RNG rng;
    for (int y = 0; y < 512; ++y)
        for (int x = 0; x < 512; ++x) {
            Point3f p(-3 + 6.f * float(x) / 512., -3 + 6.f * float(y) / 512., 0.1);
            int m = 10;
            Float sum = 0;
            for (int i = 0; i < 100; ++i)
                sum += Lewis3Stoch3(p, m, rng.Uniform<Float>());
            image.SetChannel({x, y}, 0, sum / 100);
        }

    image.Write("lewis3.exr");
}

TEST(Noise, LewisStochasticMatch) {
    RNG rng;
    for (int i = 0; i < 20; ++i) {
        Float sum[3] = {0};
        int n = 1024*1024;
        int npts = 2 + rng.Uniform<int>(8);
        auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
        Point3f p(u(), u(), u());
        Float ref = Lewis3(p, npts);
        for (int j = 0; j < n; ++j) {
            sum[0] += Lewis3Stoch1(p, npts, rng.Uniform<Float>());
            sum[1] += Lewis3Stoch2(p, npts, rng.Uniform<Float>());
            sum[2] += Lewis3Stoch3(p, npts, rng.Uniform<Float>());
        }

        EXPECT_LT(std::abs(ref - sum[0]/n), .007);
        EXPECT_LT(std::abs(ref - sum[1]/n), .07);
        EXPECT_LT(std::abs(ref - sum[2]/n), .007);
    }
}


TEST(Noise, LewisPerf) {
    RNG rng;
    std::vector<Point3f> p;
    auto u = [&]() { return -100 + 200 * rng.Uniform<Float>(); };
    for (int i = 0; i < 100; ++i)
        p.push_back(Point3f(u(), u(), u()));

    int m = 10;
    int nbench = 100000;
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += Lewis3(p[j], m);
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Normal: %.3f (sum %f)\n", elapsed, sum);
    }
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += Lewis3Stoch1(p[j], m, rng.Uniform<Float>());
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Stoch 1: %.3f (sum %f)\n", elapsed, sum);
    }
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += Lewis3Stoch2(p[j], m, rng.Uniform<Float>());
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Stoch 2: %.3f (sum %f)\n", elapsed, sum);
    }
    {
    float sum = 0;
    Timer t;
    for (int i = 0; i < nbench; ++i)
        for (int j = 0; j < p.size(); ++j)
            sum += Lewis3Stoch3(p[j], m, rng.Uniform<Float>());
    float elapsed = t.ElapsedSeconds();
    fprintf(stderr, "Stoch 3: %.3f (sum %f)\n", elapsed, sum);
    }
}
