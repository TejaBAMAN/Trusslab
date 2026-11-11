
    FROM julia:1.10

    # Set working directory
    WORKDIR /app

    # Copy Project.toml and Manifest.toml
    COPY julia_backend/Project.toml julia_backend/Manifest.toml* julia_backend/

    # Instantiate environment with low-memory mode
    RUN julia -O0 --compile=min --project=julia_backend -e 'using Pkg; Pkg.instantiate()'

    # Copy the rest of the application
    COPY . .

    # Expose port (Render sets PORT env variable)
    EXPOSE 8085

    # Start Julia server with low-memory mode (<20MB RAM, well within Render 512MB limit)
    CMD ["julia", "-O0", "--compile=min", "--project=julia_backend", "julia_backend/main.jl", "--serve"]
