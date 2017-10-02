#include <iostream>
#include <memory>

#include "solvers/md_solver.h"

/**
 * @brief Parse command line arguments
 * @param argc
 * @param argv
 * @return Path to XYZZ file
 */
std::string parse_options(int argc, char* argv[])
{
    auto getOption = [](char ** begin, char ** end, const std::string & option) -> char* {
        char ** itr = std::find(begin, end, option);

        if (itr != end && ++itr != end) {
            return *itr;
        }

        return 0;
    };

    auto optionExists = [](char** begin, char** end, const std::string& option) -> bool {
        return std::find(begin, end, option) != end;
    };

    // Search for -c flag (config file)
    if(optionExists(argv, argv + argc, "-c")) {
        std::string filename = getOption(argv, argv + argc, "-c");

        return filename;
    }
    else {
        return "";
    }
}

int main(int argc, char* argv[]) {
    auto config = parse_options(argc, argv);
    auto max_step = 1000;

    if(config.empty()) {
        std::cerr << "Usage: ./CppMD -c FILE.XYZ\n";
        return 1;
    }

    auto solver   = MDSolver::create();
         solver->read(config);

    for(std::size_t step = 0; step < max_step; ++step) {
        std::cout << " Step: " << step << std::endl;
        solver->solve();
        solver->write("step_" + std::to_string(step) + ".vtk");
    }

    return 0;
}