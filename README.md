# TrussLab: The Complete Architecture & Function Breakdown

Welcome to TrussLab! This document serves as the master blueprint for the entire `C:\random` repository. It explains every single folder, file, executable, and **every single internal function** in plain English.

---

## 1. The Core Architecture (How the App Works)

Before looking at the files, it helps to understand how the application operates:
1. **The Screen (UI):** A lightning-fast, offline web interface (HTML/JS) where you draw the bridge.
2. **The Math Brain (Julia):** A hidden engine that crunches the heavy Direct Stiffness Method (DSM) matrices.

**How they connect:**
*   **On a Windows PC:** There is **no web server**. The Python `TrussAnalyzerApp.py` script opens a native window and forces the Julia Math Brain directly into the same memory space. The Screen and the Math Brain talk to each other instantly.
*   **On an Android Phone:** Phones don't allow native memory embedding. So, the `.apk` file secretly spins up a microscopic, hidden HTTP server (`127.0.0.1:8085`) on your phone. The Screen sends a text message to this local port, the Math Brain processes it, and replies.

---

## 2. Root Folder (`C:\random`)

*   📄 **`README.md`**: The file you are reading right now!
*   🐍 **`TrussAnalyzerApp.py`**: The Windows Desktop application wrapper.
    *   `__init__()`: Bootstraps the Julia Math Brain in the background.
    *   `solve_truss()`: Intercepts the raw JSON text from the Screen and hands it directly to Julia.
    *   `main()`: Uses `pywebview` to open a blank Microsoft Edge window and injects our Screen (`index.html`) inside it, making it look like a real desktop app.
*   🐍 **`build_deb.py`**: A deployment script to package the app into a standard `.deb` Linux installer.
*   🐳 **`Dockerfile`**: A configuration file used by Docker to build the app inside an isolated virtual Linux container.

---

## 3. The `runnable_files/` Folder

This contains the actual, clickable deliverables that you distribute to users.

*   🖥️ **`TrussLabApp.exe`**: The actual native Windows Application. We used a tool called PyInstaller to freeze all the Python, HTML, and Julia code into this one single 36MB file.
*   📱 **`TrussAnalyzer.apk`**: The final Android mobile application file.
*   🐧 **`TrussAnalyzer_1.0.0_amd64.deb`**: The final Linux installer file.
*   📕 **`TrussLab_Architecture_Detailed.pdf`**: A highly detailed, formatted PDF document containing a highly technical breakdown for software engineers.

---

## 4. The `apk/` Folder (Android Builder)

Contains the raw tools required to compile the Android app. 

*   🐍 **`build_apk.py`**: The master automation script. Running this grabs the UI, triggers the `aapt.exe` and `zipalign.exe` tools, and spits out the final `.apk` file.
*   ⚙️ **`AndroidManifest.xml`**: Tells the Android phone what the app is called, its icon, and requests permissions (like internet access).
*   🛠️ **`aapt.exe`**: (Android Asset Packaging Tool). Takes raw images and XML files and compiles them into a dense binary format.
*   🛠️ **`zipalign.exe`**: Optimizes the memory alignment of the final `.apk` file so the app uses less RAM on the phone.
*   📚 **`android.jar`**: The core Android Operating System library required for compilation.
*   📁 **`res/mipmap-*/`**: Folders holding the `ic_launcher.png` app icons for different phone screen resolutions (hdpi, xhdpi, xxhdpi).

---

## 5. The `apk/assets/www/` Folder (The Screen)

*   🌐 **`index.html`**: The single, monolithic UI file. It contains all HTML, CSS, and Javascript.
    *   `newProject()` / `showHome()`: Cleans the UI slate to start a new bridge or return to the main menu.
    *   `saveProject()` / `editProject()`: Manages saving drafts to the browser's hidden memory and loading them back.
    *   `go(step)`: The traffic cop. It moves you sequentially from Step 1 through Step 7.
    *   `validateStep(s)`: The security guard. Blocks you from moving to the next step if you draw an impossible triangle or forget a material property.
    *   `buildTrianglesUI()`: Dynamically creates the HTML input boxes based on how many triangles you asked for.
    *   `generateGeometry()`: Takes your lengths (e.g., "5 meters") and uses the Law of Cosines to calculate exact X and Y coordinates on a graph.
    *   `drawGeometry()` / `drawTruss()` / `drawLoad()`: The artists. They physically paint the blue structural lines, nodes, and red arrows onto the HTML5 `<canvas>`.
    *   `solve()`: Packages up your final drawing into a JSON text string and fires it off to the Math Brain.

---

## 6. The `julia_backend/` Folder (The Math Brain)

The scientific physics engine, split into specific modules.

🧠 **`solver.jl`**
*   `solve_truss_fem()`: The ultimate manager function. It receives the bridge from the UI, delegates the math to the functions below, and returns the final answers.

🏗️ **`assembly.jl`**
*   `assemble_global_stiffness()`: Looks at every single steel beam and mathematically glues them together into one giant interconnected puzzle (the Global Stiffness Matrix [K]).
*   `assemble_global_load_vector()`: Gathers all the weights/forces you applied and puts them into a single mathematical list [F].

📐 **`stiffness.jl`**
*   `calculate_local_stiffness()`: Calculates exactly how stiff a single steel beam is natively (EA/L).
*   `calculate_element_stiffness()`: Rotates that stiffness based on the beam's angle (using sine/cosine).

⚓ **`boundary_conditions.jl`**
*   `get_node_dofs()`: Figures out the exact X and Y movement directions for a specific joint.
*   `identify_dofs()`: Finds which nodes are bolted to the ground (Pins/Rollers) and locks those specific X/Y directions so the bridge doesn't float away.

🗣️ **`communication.jl`**
*   `handle_request_json()` & `parse_and_process()`: Translates the raw text JSON received from the Screen into strict Julia variables.
*   `dict_to_triangles()`, `dict_to_supports()`, `dict_to_loads()`: Maps the parsed text into physical structural objects.

🧩 **`elements.jl`**
*   `generate_unique_elements()`: Algorithmically scans adjacent triangles and deletes duplicate edges, finding the true number of unique steel beams.

🗺️ **`geometry.jl`**
*   `calculate_triangle_angles()` & `distance_2d()`: Performs heavy trigonometry to double-check angles and distances.
*   `calculate_edge_orientation()`: Figures out the physical slope/angle of a beam.
*   `triangle_centroid()`, `signed_triangle_area()`, `point_in_triangle()`: Advanced geometric mapping functions.

🗺️ **`triangles.jl`**
*   `find_shared_edge_info()` & `reconstruct_geometry_from_triangles()`: Backend verification of the frontend's coordinate mapping to ensure no data tampering occurred.

👷 **`validation.jl`**
*   `validate_triangle_lengths()`: Proves that your triangles are physically possible (Side A + Side B > Side C).
*   `check_zero_length()`: Makes sure you didn't accidentally create a steel beam with a length of 0.
*   `validate_connectivity()` & `check_duplicate_nodes()`: Ensures all parts of the bridge are actually connected to each other.

⚖️ **`stability.jl`**
*   `check_structural_stability()`: Runs kinematic determinacy checks (b + r >= 2j) to ensure the structure isn't a free-moving, collapsing mechanism before attempting to solve.

📊 **`results.jl`**
*   `compute_deformed_coordinates()`: Calculates exactly how far each node sags under the weight so the Screen can draw the bent bridge.
*   `generate_summary()`: Calculates the internal stress (tension/compression) inside each beam and the reaction forces hitting the anchors.

🚀 **`main.jl` & `Project.toml`**
*   `start_server()` & `run_cli()`: The standalone entry points to boot up the invisible Mobile HTTP server or run the Math Engine from a terminal.
*   `Project.toml`: Julia's package manager file. Ensures dependencies like `JSON3` are downloaded.

🧪 **`tests/runtests.jl`**
*   Contains automated testing blocks to ensure if we change a math equation, it doesn't accidentally break the physics engine.
