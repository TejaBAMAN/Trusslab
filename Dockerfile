FROM julia:1.10

# Set working directory
WORKDIR /app

# Copy Project.toml and Manifest.toml for package installation
COPY julia_backend/Project.toml julia_backend/Manifest.toml* julia_backend/

# Instantiate environment with explicit project flag
RUN julia --project=julia_backend -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy the rest of the application
COPY . .

# Expose port (Render sets PORT env variable)
EXPOSE 8085

# Start the Julia server with explicit project flag
CMD ["julia", "--project=julia_backend", "julia_backend/main.jl", "--serve"]
