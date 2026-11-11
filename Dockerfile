FROM julia:1.10

# Set working directory
WORKDIR /app

# Copy Project.toml first for caching
COPY julia_backend/Project.toml julia_backend/

# Instantiate environment
RUN cd julia_backend && julia -e 'using Pkg; Pkg.instantiate()'

# Copy the rest of the application
COPY . .

# Expose port (Render sets PORT env variable)
EXPOSE 8085

# Start the Julia server
CMD ["julia", "julia_backend/main.jl", "--serve"]
