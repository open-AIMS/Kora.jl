using BenchmarkTools
using Kora
using Random

const SUITE = BenchmarkGroup()

# --- Initialization benchmarks ---
SUITE["initialization"] = BenchmarkGroup()

SUITE["initialization"]["initialize_reef"] = @benchmarkable begin
    Kora.initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0)
end

SUITE["initialization"]["initialize_population"] = @benchmarkable begin
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0)
end

# --- Core simulation benchmarks ---
SUITE["simulation"] = BenchmarkGroup()

SUITE["simulation"]["run_model_small"] = @benchmarkable begin
    Kora.reset!(reef)
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
    Kora.run_model!(reef, env; rng=Xoshiro(1))
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=20, n_locs=3, density=20, area=90.0)
    env = Kora.generate_example_environment(20, 3; rng=Xoshiro(42))
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
end

SUITE["simulation"]["run_model_medium"] = @benchmarkable begin
    Kora.reset!(reef)
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
    Kora.run_model!(reef, env; rng=Xoshiro(1))
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0)
    env = Kora.generate_example_environment(50, 20; rng=Xoshiro(42))
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
end

SUITE["simulation"]["run_model_large"] = @benchmarkable begin
    Kora.reset!(reef)
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
    Kora.run_model!(reef, env; rng=Xoshiro(1))
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=75, n_locs=100, density=20, area=90.0)
    env = Kora.generate_example_environment(75, 100; rng=Xoshiro(42))
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
end

# --- Cover calculation benchmarks ---
SUITE["cover"] = BenchmarkGroup()

SUITE["cover"]["coral_cover_total"] = @benchmarkable begin
    Kora.coral_cover(reef)
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0)
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
    env = Kora.generate_example_environment(50, 20; rng=Xoshiro(42))
    Kora.run_model!(reef, env; rng=Xoshiro(1))
end

SUITE["cover"]["group_cover"] = @benchmarkable begin
    Kora.group_cover(reef)
end setup = begin
    reef = Kora.initialize_reef(; n_timesteps=50, n_locs=20, density=20, area=90.0)
    Kora.initialize_coral_population!(reef; rng=Xoshiro(1))
    env = Kora.generate_example_environment(50, 20; rng=Xoshiro(42))
    Kora.run_model!(reef, env; rng=Xoshiro(1))
end

# --- Component benchmarks ---
SUITE["components"] = BenchmarkGroup()

SUITE["components"]["cover_cm_to_m2"] = @benchmarkable begin
    Kora.cover_cm_to_m2.(diams)
end setup = begin
    diams = rand(Float32, 1000) .* 100.0f0
end

SUITE["components"]["bleaching_susceptibility"] = @benchmarkable begin
    Kora.bleaching_susceptibility.(diams)
end setup = begin
    diams = rand(Float32, 1000) .* 200.0f0
end

SUITE["components"]["generate_environment"] = @benchmarkable begin
    Kora.generate_example_environment(75, 100; rng=Xoshiro(42))
end
