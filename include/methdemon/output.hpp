#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include "tumour.hpp"
#include <fstream>
#include <iostream>
#include <string>

class FileOutput {
private:
    std::ofstream file;
public:
    // Constructor and destructor
    FileOutput(const std::string& path) { file.open(path, std::ofstream::out); }
    ~FileOutput() { file.close(); }
    // Write to file
    void writeDemesFile(Tumour& tumour, float gensElapsed);
    void writeDemesHeader();
};

#endif // OUTPUT_HPP