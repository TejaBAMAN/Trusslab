# julia_backend/main.jl
import Pkg
Pkg.activate(@__DIR__)
# Main entry point for Julia Truss Engineering and Analysis Backend

include("validation.jl")
include("geometry.jl")
include("triangles.jl")
include("elements.jl")
include("stiffness.jl")
include("boundary_conditions.jl")
include("assembly.jl")
include("stability.jl")
include("results.jl")
include("solver.jl")
include("communication.jl")

using .Validation
using .Geometry
using .Triangles
using .Elements
using .Stiffness
using .BoundaryConditions
using .Assembly
using .Stability
using .Results
using .Solver
using .Communication

using HTTP
using JSON
using Printf

function start_server(port::Int=8085)
    println("="^60)
    println("JULIA 2D TRUSS ENGINEERING & FEM ANALYSIS ENGINE")
    println("Engine version: 1.0.0")
    println("Listening on http://127.0.0.1:$(port)")
    println("="^60)
    
    router = HTTP.Router()
    
    # Health check endpoint
    HTTP.register!(router, "GET", "/api/ping", function(req::HTTP.Request)
        return HTTP.Response(200, [
            "Content-Type" => "application/json",
            "Access-Control-Allow-Origin" => "*",
            "Access-Control-Allow-Methods" => "GET, OPTIONS",
            "Access-Control-Allow-Headers" => "Content-Type"
        ],
            JSON.json(Dict("status" => "online", "engine" => "Julia Truss FEM Engine", "version" => "1.0.0"))
        )
    end)

    HTTP.register!(router, "OPTIONS", "/api/ping", function(req::HTTP.Request)
        return HTTP.Response(200, [
            "Access-Control-Allow-Origin" => "*",
            "Access-Control-Allow-Methods" => "GET, OPTIONS",
            "Access-Control-Allow-Headers" => "Content-Type"
        ], "")
    end)
    
    # Primary API endpoint for analysis and geometry reconstruction
    HTTP.register!(router, "POST", "/api/truss", function(req::HTTP.Request)
        try
            body_str = String(req.body)
            resp_str = handle_request_json(body_str)
            return HTTP.Response(200, [
                "Content-Type" => "application/json",
                "Access-Control-Allow-Origin" => "*",
                "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                "Access-Control-Allow-Headers" => "Content-Type"
            ], resp_str)
        catch err
            err_dict = Dict(
                "success" => false,
                "error" => "Internal Julia Server Error: $(err)",
                "stacktrace" => sprint(showerror, err, catch_backtrace())
            )
            return HTTP.Response(500, ["Content-Type" => "application/json"], JSON.json(err_dict))
        end
    end)
    
    # Handle CORS preflight
    HTTP.register!(router, "OPTIONS", "/api/truss", function(req::HTTP.Request)
        return HTTP.Response(200, [
            "Access-Control-Allow-Origin" => "*",
            "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
            "Access-Control-Allow-Headers" => "Content-Type"
        ], "")
    end)
    
    # Serve Web UI (Mobile/Tablet/PC)
    www_dir = normpath(joinpath(@__DIR__, "..", "www"))
    if !isfile(joinpath(www_dir, "index.html"))
        www_dir = normpath(joinpath(@__DIR__, "..", "apk", "assets", "www"))
    end
    
    HTTP.register!(router, "GET", "/", function(req::HTTP.Request)
        html_path = joinpath(www_dir, "index.html")
        if isfile(html_path)
            return HTTP.Response(200, ["Content-Type" => "text/html"], read(html_path))
        else
            return HTTP.Response(404, "Frontend not found at $html_path")
        end
    end)

    HTTP.register!(router, "GET", "/app.js", function(req::HTTP.Request)
        js_path = joinpath(www_dir, "app.js")
        if isfile(js_path)
            return HTTP.Response(200, ["Content-Type" => "application/javascript"], read(js_path))
        else
            return HTTP.Response(404, "app.js not found")
        end
    end)

    HTTP.register!(router, "GET", "/style.css", function(req::HTTP.Request)
        css_path = joinpath(www_dir, "style.css")
        if isfile(css_path)
            return HTTP.Response(200, ["Content-Type" => "text/css"], read(css_path))
        else
            return HTTP.Response(404, "style.css not found")
        end
    end)
    
    HTTP.serve(router, "0.0.0.0", port)
end

function run_cli(args::Vector{String})
    if "--help" in args || "-h" in args
        println("""
        Usage: julia main.jl [options]
        
        Options:
          --serve [port]              Run HTTP server (default port: 8085)
          --input <in.json> --output <out.json>
                                      Process single JSON request file
          --test                      Run backend test suite
          -h, --help                  Show this help message
        """)
        return
    end
    
    if "--test" in args
        include(joinpath(@__DIR__, "tests", "runtests.jl"))
        return
    end
    
    # Check for --serve
    serve_idx = findfirst(x -> x == "--serve", args)
    if serve_idx !== nothing
        port = 8085
        if haskey(ENV, "PORT")
            port = parse(Int, ENV["PORT"])
        end
        if serve_idx < length(args)
            p_arg = tryparse(Int, args[serve_idx + 1])
            if p_arg !== nothing
                port = p_arg
            end
        end
        start_server(port)
        return
    end
    
    # Check for --input and --output
    in_idx = findfirst(x -> x == "--input", args)
    out_idx = findfirst(x -> x == "--output", args)
    if in_idx !== nothing && in_idx < length(args)
        in_file = args[in_idx + 1]
        out_file = (out_idx !== nothing && out_idx < length(args)) ? args[out_idx + 1] : "output.json"
        
        println("Processing JSON request from: $in_file")
        input_str = read(in_file, String)
        out_str = handle_request_json(input_str)
        write(out_file, out_str)
        println("Results written to: $out_file")
        return
    end
    
    # Default behavior if no arguments: start server on 8085
    start_server(8085)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_cli(ARGS)
end
