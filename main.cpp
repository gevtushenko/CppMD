#include <iostream>
#include <iomanip>
#include <memory>
#include <chrono>

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

void print_step_info(const std::size_t& step, const double& time, const double& dt, const double& real_time)
{
    std::cout << "  \033[1;35mStep\033[0m:"  << std::setw(4)  << step
              << "; \033[1;35mTime\033[0m: " << std::setw(10) << time
              << "; \033[1;35mDT\033[0m: "   << std::setw(10) << dt
              << "; \033[1;35mIn\033[0m: "   << std::setw(10) << real_time << "s\n";
}

int main(int argc, char* argv[]) {
    auto config = parse_options(argc, argv);
    auto max_step = 1000;
    auto write_every = 20;

    if(config.empty()) {
        std::cerr << "Usage: ./CppMD -c FILE.XYZ\n";
        return 1;
    }

    auto solver   = MDSolver::create();
         solver->read(config);

    double time = 0.0;

    for(std::size_t step = 0; step < max_step; ++step) {
        auto step_begin = std::chrono::steady_clock::now();

        solver->solve(); time += solver->dt();

        auto   step_end = std::chrono::steady_clock::now();
        double step_time = std::chrono::duration<double>(step_end - step_begin).count();
        print_step_info(step, time, solver->dt(), step_time );

        if(step & write_every == 0) {
            solver->write("step_" + std::to_string(step) + ".vtk");
        }
    }

    return 0;
}