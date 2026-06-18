#include "DataExporter.h"

#include <fstream>
#include <iostream>

namespace DataExporter {
    bool exportToObj(const MACGrid& grid, const std::string& filepath) {
        std::ofstream outfile(filepath);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            return false;
        }

        const auto& particles = grid.particles();
        outfile << "# Particle data exported from FluidSim\n";
        outfile << "# Total particles: " << particles.size() << "\n";
        for (const auto& particle : particles) {
            outfile << "v " << particle.position.x << " " << particle.position.y << " " << particle.position.z << "\n";
        }

        return true;
    }

    bool exportMeshToObj(const TriangleMesh& mesh, const std::string& filepath) {
        std::ofstream outfile(filepath);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            return false;
        }

        outfile << "# Mesh data exported from FluidSim\n";
        outfile << "# Vertices: " << mesh.vertices.size() << "\n";
        outfile << "# Faces: " << mesh.faces.size() / 3 << "\n";

        for (const auto& v : mesh.vertices) {
            outfile << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }

        for (const auto& n : mesh.normals) {
            outfile << "vn " << n.x << " " << n.y << " " << n.z << "\n";
        }

        for (size_t i = 0; i < mesh.faces.size(); i += 3) {
            outfile << "f " << mesh.faces[i] + 1 << " "
                    << mesh.faces[i + 1] + 1 << " "
                    << mesh.faces[i + 2] + 1 << "\n";
        }

        return true;
    }
}
