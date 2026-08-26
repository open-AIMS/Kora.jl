# build/kora_server_main.jl -- juliac --output-exe entry point for kora-server
#
# Untrimmed juliac build (matches ReefGuideWorker.jl/build/worker_main.jl's
# shape -- no --trim, AOT-compiles what's statically inferable and keeps the
# dynamic compiler embedded for the rest). Oxygen's macro-based routing and
# HTTP.serve's internals rely on this: a Phase 2 spike (a minimal Oxygen
# "hello world" compiled the same way) confirmed an untrimmed
# `juliac --output-exe --bundle` build starts and serves real HTTP traffic.
#
# Compiled against build/server/Project.toml (Kora as a path dep, Oxygen/HTTP
# as real deps) rather than Kora.jl's own top-level Project.toml, so Oxygen
# support stays behind Kora's weakdep extension mechanism (ext/OxygenExt.jl)
# for every other consumer (kora-worker, the native desktop bridge) -- see
# .claude/plans/web-app/kora-web-service.md Component 1.
#
# All server behavior (routes, session registry, idle timeout, env-var
# config) lives in ext/OxygenExt.jl / Kora.start_server -- this file is a
# thin juliac entry point only.

using Kora
using Oxygen  # loads OxygenExt, giving Kora.start_server a method

function main(ARGS::Vector{String})::Cint
    Kora.start_server()
    return Cint(0)
end

Base.@main
