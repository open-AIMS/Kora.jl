using Test
using Kora
using Random
using Statistics

const RNG = Xoshiro(42)
const MODELS_AVAILABLE = !isnothing(Kora.growth_models) && !isnothing(Kora.survival_models)

@testset "Unit Tests" begin
    @testset "Statistical Functions" begin
        @testset "rational_erf" begin
            @test Kora.rational_erf(0.0f0) == 0.0f0
            @test Kora.rational_erf(1.0f0) > 0.0f0
            @test Kora.rational_erf(-1.0f0) < 0.0f0
            # erf is odd: erf(-x) == -erf(x)
            @test Kora.rational_erf(-0.5f0) ≈ -Kora.rational_erf(0.5f0) atol = 1e-5
            # erf(large) -> 1
            @test Kora.rational_erf(4.0f0) ≈ 1.0f0 atol = 1e-4
        end

        @testset "truncated_normal_mean" begin
            # Standard normal truncated symmetrically should have mean ~0
            m = Kora.truncated_standard_normal_mean(-2.0f0, 2.0f0)
            @test abs(m) < 0.01f0

            # Truncated to positive half should have positive mean
            m_pos = Kora.truncated_standard_normal_mean(0.0f0, 3.0f0)
            @test m_pos > 0.0f0

            # General truncated normal
            m_gen = Kora.truncated_normal_mean(5.0f0, 1.0f0, 3.0f0, 7.0f0)
            @test 3.0f0 < m_gen < 7.0f0
            @test abs(m_gen - 5.0f0) < 0.5f0  # close to the untruncated mean
        end

        @testset "truncated_normal_cdf" begin
            # CDF at lower bound should be 0
            @test Kora.truncated_normal_cdf(0.0f0, 0.0f0, 1.0f0, 0.0f0, 2.0f0) ≈ 0.0f0 atol =
                1e-5
            # CDF at upper bound should be 1
            @test Kora.truncated_normal_cdf(2.0f0, 0.0f0, 1.0f0, 0.0f0, 2.0f0) ≈ 1.0f0 atol =
                1e-3
            # CDF should be monotonically increasing
            c1 = Kora.truncated_normal_cdf(0.5f0, 0.0f0, 1.0f0, -2.0f0, 2.0f0)
            c2 = Kora.truncated_normal_cdf(1.0f0, 0.0f0, 1.0f0, -2.0f0, 2.0f0)
            @test c2 > c1
        end
    end

    @testset "Metrics" begin
        y = Float32[1.0, 2.0, 3.0, 4.0, 5.0]
        y_hat = Float32[1.1, 2.0, 2.9, 4.1, 5.0]

        @test Kora.RMSE(y_hat, y) < 0.1
        @test Kora.R2(y_hat, y) > 0.99
        @test Kora.pearson(y_hat, y) > 0.99
        @test Kora.spearman(y_hat, y) > 0.99
        @test Kora.kendall(y_hat, y) > 0.9

        # Perfect predictions
        @test Kora.RMSE(y, y) == 0.0
        @test Kora.R2(y, y) ≈ 1.0
    end

    @testset "Size Classes and Cover" begin
        @testset "bin_edges" begin
            edges = Kora.bin_edges()
            @test size(edges, 1) == length(Kora.TARGET_GROUPS)
            @test size(edges, 2) > 0
            # Bin edges should be non-negative and increasing per group
            for g in 1:size(edges, 1)
                row = edges[g, :]
                valid = row[row .> 0]
                @test issorted(valid)
            end
        end

        @testset "cover_cm_to_m2" begin
            # A coral with diameter 0 has 0 cover
            @test Kora.cover_cm_to_m2(0.0f0) == 0.0f0
            # Cover should be positive for positive diameter
            @test Kora.cover_cm_to_m2(10.0f0) > 0.0f0
            # Cover should scale with diameter squared (circle area)
            c1 = Kora.cover_cm_to_m2(10.0f0)
            c2 = Kora.cover_cm_to_m2(20.0f0)
            @test c2 / c1 ≈ 4.0f0 atol = 0.01
        end

        @testset "thresholds" begin
            sus = Kora.susceptibility_size_thresholds()
            mat = Kora.mature_size_thresholds()
            @test length(sus) == length(Kora.TARGET_GROUPS)
            @test length(mat) == length(Kora.TARGET_GROUPS)
            @test all(sus .> 0)
            @test all(mat .> 0)
        end
    end

    @testset "Coral Mortality" begin
        @testset "bleaching_susceptibility" begin
            # Small corals are more susceptible (size refuge effect)
            s_small = Kora.bleaching_susceptibility(5.0f0)
            s_large = Kora.bleaching_susceptibility(200.0f0)
            @test s_small > s_large
            # Output in [0, 1]
            @test 0.0f0 <= s_small <= 1.0f0
            @test 0.0f0 <= s_large <= 1.0f0
        end

        @testset "depth_coefficient" begin
            # Shallow depth -> higher bleaching
            shallow = Kora.depth_coefficient(2.0f0)
            deep = Kora.depth_coefficient(15.0f0)
            @test shallow > deep
            @test shallow > 0.0f0
            @test deep > 0.0f0
        end
    end

    @testset "Recruitment" begin
        @testset "larval_production" begin
            # Larger corals produce more larvae
            lp_small = Kora.larval_production(5.0f0, 1)
            lp_large = Kora.larval_production(50.0f0, 1)
            @test lp_large > lp_small
            @test lp_small >= 0.0f0
        end

        @testset "breeders equation" begin
            # No selection differential -> no change
            @test Kora.breeders(5.0f0, 5.0f0, 0.5f0) ≈ 5.0f0 atol = 1e-5
            # Positive selection: offspring mean shifts toward selected mean
            result = Kora.breeders(5.0f0, 7.0f0, 0.5f0)
            @test result > 5.0f0
            @test result < 7.0f0
        end
    end
end

@testset "Fitted Models" begin
    @testset "Growth Model" begin
        gm = Kora.growth_models
        if isnothing(gm)
            @warn "Skipping growth model tests: pre-fitted models not available"
        else
            @test length(gm) == length(Kora.TARGET_GROUPS)
            for g in 1:length(Kora.TARGET_GROUPS)
                growth = gm[g](10.0f0)
                @test isfinite(growth)
            end
        end
    end

    @testset "Survival Model" begin
        sm = Kora.survival_models
        if isnothing(sm)
            @warn "Skipping survival model tests: pre-fitted models not available"
        else
            @test length(sm) == length(Kora.TARGET_GROUPS)
            for g in 1:length(Kora.TARGET_GROUPS)
                surv = sm[g](10.0f0)
                @test isfinite(surv)
                @test 0.0f0 <= surv <= 1.0f0
            end
        end
    end
end

@testset "ReefState and Environment" begin
    @testset "ReefState Initialization" begin
        if isnothing(Kora.growth_models) || isnothing(Kora.survival_models)
            @warn "Skipping ReefState tests: pre-fitted models not available"
        else
            rng = Xoshiro(42)
            n_ts = 10
            n_locs = 5

            reef = Kora.initialize_reef(;
                n_timesteps=n_ts,
                n_locs=n_locs,
                density=20,
                area=90.0
            )

            @test Kora.n_timesteps(reef) == n_ts
            @test Kora.n_locations(reef) == n_locs
            @test Kora.n_groups(reef) == length(Kora.TARGET_GROUPS)

            # Before initialization, populations should be empty
            @test Kora.total_population(reef, 1, 1) == 0

            # Initialize populations
            Kora.initialize_coral_population!(reef; rng=rng)

            # After initialization, each location should have corals
            for loc in 1:n_locs
                @test Kora.total_wild(reef, 1, loc) > 0
            end

            # Deployed population starts at zero
            for loc in 1:n_locs
                @test Kora.total_deployed(reef, 1, loc) == 0
            end
        end
    end

    @testset "Environment Generation" begin
        rng = Xoshiro(123)
        n_years = 20
        n_locs = 5

        env = Kora.generate_example_environment(n_years, n_locs; rng=rng)

        # Should have correct dimensions
        @test size(env, 1) == n_years  # time axis
        @test size(env, 2) == n_locs   # location axis

        # DHW values should be non-negative
        @test all(env .>= 0.0f0)
    end
end

@testset "Integration" begin
    @testset "Full Simulation Run" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            rng = Xoshiro(99)
            n_ts = 20
            n_locs = 3

            reef = Kora.initialize_reef(;
                n_timesteps=n_ts,
                n_locs=n_locs,
                density=20,
                area=90.0
            )
            Kora.initialize_coral_population!(reef; rng=rng)

            env = Kora.generate_example_environment(n_ts, n_locs; rng=rng)

            # Run simulation and log elapsed time
            elapsed = @elapsed Kora.run_model!(reef, env; rng=rng)
            @info "Full simulation run" elapsed_seconds = round(elapsed; digits = 3) n_timesteps = n_ts n_locations = n_locs

            # After simulation, cover should be computed
            cover = Kora.coral_cover(reef)
            @test length(cover) == n_ts
            @test all(isfinite.(cover))
            @test all(cover .>= 0.0f0)

            # Initial cover should be positive (we initialized populations)
            @test cover[1] > 0.0f0

            # Group cover should have correct dimensions
            gc = Kora.group_cover(reef)
            @test size(gc, 1) == n_ts
            @test size(gc, 2) == length(Kora.TARGET_GROUPS)

            # Location-level cover
            loc_cover = Kora.coral_cover(reef, 1)
            @test length(loc_cover) == n_locs
            @test all(loc_cover .>= 0.0f0)
        end
    end

    @testset "Reproducibility" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            function run_with_seed(seed)
                rng = Xoshiro(seed)
                reef = Kora.initialize_reef(; n_timesteps=15, n_locs=2, density=20, area=90.0)
                Kora.initialize_coral_population!(reef; rng=rng)
                env = Kora.generate_example_environment(15, 2; rng=Xoshiro(seed))
                Kora.run_model!(reef, env; rng=Xoshiro(seed))
                return Kora.coral_cover(reef)
            end

            cover_a = run_with_seed(777)
            cover_b = run_with_seed(777)
            @test cover_a == cover_b

            # Different seed should (almost certainly) produce different results
            cover_c = run_with_seed(888)
            @test cover_a != cover_c
        end
    end

    @testset "Reset" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            rng = Xoshiro(55)
            reef = Kora.initialize_reef(; n_timesteps=10, n_locs=2, density=20, area=90.0)
            Kora.initialize_coral_population!(reef; rng=rng)
            env = Kora.generate_example_environment(10, 2; rng=rng)

            Kora.run_model!(reef, env; rng=rng)

            # Reset should clear simulation results
            Kora.reset!(reef)

            # After reset, only timestep 1 should have population
            @test Kora.total_population(reef, 1, 1) > 0
            # Later timesteps should be empty
            @test Kora.total_population(reef, 5, 1) == 0
        end
    end

    @testset "Space Constraint" begin
        # space_constraint should limit growth as cover approaches capacity
        low = Kora.space_constraint(0.1f0, 1.0f0)
        high = Kora.space_constraint(0.9f0, 1.0f0)
        @test low > high  # More space available at low cover
        @test low <= 1.0f0
        @test high >= 0.0f0
    end
end

@testset "Model I/O" begin
    @testset "Round-trip: growth model" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            gm = Kora.growth_models
            path = tempname() * ".json"
            try
                Kora.save_models(gm, path)
                gm2 = Kora.load_models(path)
                @test gm2 isa Kora.PolyGrowthModel
                @test length(gm2) == length(gm)
                for g in 1:length(gm)
                    @test gm[g](10.0f0) == gm2[g](10.0f0)
                    @test gm[g](50.0f0) == gm2[g](50.0f0)
                end
            finally
                isfile(path) && rm(path)
            end
        end
    end

    @testset "Round-trip: survival model" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            sm = Kora.survival_models
            path = tempname() * ".json"
            try
                Kora.save_models(sm, path)
                sm2 = Kora.load_models(path)
                @test sm2 isa Kora.PolySurvivalModel
                @test length(sm2) == length(sm)
                for g in 1:length(sm)
                    @test sm[g](10.0f0) == sm2[g](10.0f0)
                    @test sm[g](50.0f0) == sm2[g](50.0f0)
                end
            finally
                isfile(path) && rm(path)
            end
        end
    end

    @testset "load_models: unsupported format_version" begin
        path = tempname() * ".json"
        try
            write(path, """
            {
              "format_version": 999,
              "model_kind": "growth",
              "fitted_at": "2024-01-01T00:00:00",
              "models": [],
              "performance": {"train": {}, "test": {}}
            }
            """)
            @test_throws Exception Kora.load_models(path)
        finally
            isfile(path) && rm(path)
        end
    end

    @testset "load_models: unknown model_kind" begin
        path = tempname() * ".json"
        try
            write(path, """
            {
              "format_version": 1,
              "model_kind": "unknown_kind",
              "fitted_at": "2024-01-01T00:00:00",
              "models": [],
              "performance": {"train": {}, "test": {}}
            }
            """)
            @test_throws Exception Kora.load_models(path)
        finally
            isfile(path) && rm(path)
        end
    end

    @testset "load_models: missing required top-level fields" begin
        path = tempname() * ".json"
        try
            write(path, """{"format_version": 1}""")
            @test_throws Exception Kora.load_models(path)
        finally
            isfile(path) && rm(path)
        end
    end

    @testset "load_models: unknown model type tag" begin
        path = tempname() * ".json"
        try
            write(path, """
            {
              "format_version": 1,
              "model_kind": "growth",
              "fitted_at": "2024-01-01T00:00:00",
              "models": [
                {"type": "UnknownModelType", "name": "test",
                 "dtype": "float32", "min_x": 1.0, "min_y": 0.5,
                 "max_x": 100.0, "max_y": 10.0, "poly_coeffs": [1.0, 0.5]}
              ],
              "performance": {"train": {}, "test": {}}
            }
            """)
            @test_throws Exception Kora.load_models(path)
        finally
            isfile(path) && rm(path)
        end
    end

    @testset "check_model_pair_skew: matching timestamps" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            gm = Kora.growth_models
            sm = Kora.survival_models
            g_path = tempname() * ".json"
            s_path = tempname() * ".json"
            try
                Kora.save_models(gm, g_path)
                Kora.save_models(sm, s_path)
                # Same timestamp -> no warning expected (just shouldn't throw)
                @test_nowarn Kora.check_model_pair_skew(g_path, g_path)
            finally
                isfile(g_path) && rm(g_path)
                isfile(s_path) && rm(s_path)
            end
        end
    end

    @testset "check_model_pair_skew: large skew warns" begin
        g_path = tempname() * ".json"
        s_path = tempname() * ".json"
        try
            write(g_path, """{"fitted_at": "2024-01-01T00:00:00"}""")
            write(s_path, """{"fitted_at": "2025-06-01T00:00:00"}""")
            @test_logs (:warn,) Kora.check_model_pair_skew(
                g_path, s_path; threshold_seconds=86400
            )
        finally
            isfile(g_path) && rm(g_path)
            isfile(s_path) && rm(s_path)
        end
    end
end

@testset "Regression" begin
    @testset "Cover Bounds" begin
        if !MODELS_AVAILABLE
            @warn "Skipping: pre-fitted models not available"
        else
            rng = Xoshiro(2024)
            reef = Kora.initialize_reef(; n_timesteps=50, n_locs=5, density=20, area=90.0)
            Kora.initialize_coral_population!(reef; rng=rng)
            env = Kora.generate_example_environment(50, 5; rng=rng)
            Kora.run_model!(reef, env; rng=rng)

            for ts in 1:50, loc in 1:5
                c = Kora.coral_cover(reef, ts, loc)
                @test c >= 0.0f0
                # Cover (m²) must not exceed carrying capacity (m²)
                @test c <= reef.carrying_capacity[loc]
            end
        end
    end
end
