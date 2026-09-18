
#pragma once

#include <string>

#include "ordering.hpp"


struct FactorResult {
    bool ok = false;          
    bool threw = false;        
    std::string matrix;        
    std::string error;        
    std::string summary;       
    std::string report;        
    std::string suggested_file_name;
};

FactorResult factorize_matrix(const std::string& path, Ordering ordering);
